# Usage
#
# This script is multi-threaded, start with
#     julia --project=. -t 16 --gcthreads=12,1
#
# include("scripts/validate-with-mlj.jl")
#
# # Repeat each fold 1 time:
# @time res, fails, changes = validate_all(; rng_seed=42, repeats=1)
#
# # Repeat each fold 100 times:
# @time res, fails, changes = validate_all(; rng_seed=42, repeats=100)

using DataFrames: DataFrame, nrow
using Dates: Date
using StatsBase: mean
using StableRNGs: StableRNG
import TOML
import MLJ
import MLJBase
import MLJModelInterface
import ScientificTypesBase
const MMI = MLJModelInterface

import AquaCrop
import CropGrowthTutorial

include("willmott.jl")
using .Willmott: willmott_d


const DEFAULT_MEASURES = [MLJ.rms, MLJ.mae, MLJ.mape, MLJ.rmsp, willmott_d]
const DEFAULT_MAXITER = 100

# MLJ wrapper type with no parameters tunable via MLJ, needs custom
# fitting function
MMI.@mlj_model mutable struct AquaCropFit <: MLJModelInterface.Deterministic
    rng_seed = rand(Int)
    maxiter = DEFAULT_MAXITER
end
MMI.input_scitype(::Type{<:AquaCropFit}) = ScientificTypesBase.Unknown
MMI.target_scitype(::Type{<:AquaCropFit}) = ScientificTypesBase.Unknown

# Inputs needed for model fitting for AquaCropFit
@kwdef struct CropInputs
    station_name::String
    crop_type::String
    aquacrop_name::String
    soil_type::String
    planting_density::Int
    weather_df::DataFrame
    phenology_raw_df::DataFrame
    years::Vector{Int}
end
Base.length(X::CropInputs) = length(X.years)
MLJBase.nrows(X::CropInputs) = length(X.years)
function MLJBase.selectrows(X::CropInputs, rows)
    years = X.years[rows]
    # We add one extra year, so that the phenology data for the
    # harvest date (which is in the following year) is included
    phenology_raw_df = filter(row -> row.referenzjahr in [years..., years[end]+1], X.phenology_raw_df)
    return CropInputs(;
        station_name = X.station_name,
        crop_type = X.crop_type,
        aquacrop_name = X.aquacrop_name,
        soil_type = X.soil_type,
        planting_density = X.planting_density,
        weather_df = X.weather_df,
        phenology_raw_df,
        years,
    )
end

# helper function
function get_fit_type_and_phases(crop_type::AbstractString)
    if startswith(crop_type, "silage")
        # silage crops, fit to biomass
        fit_type = :Biomass
        sowing_phase = 10  # numerical id of the sowing phase
        harvest_phase = 39  # numerical id of the harvest phase
    else
        # non-silage crops, fit to yield
        fit_type = :Yield
        sowing_phase = 10  # numerical id of the sowing phase
        harvest_phase = 24  # numerical id of the harvest phase
    end
    return fit_type, sowing_phase, harvest_phase
end

function MMI.fit(model::AquaCropFit, verbosity::Int, X::CropInputs, y::AbstractVector)
    # returns (fitresult::Dict, cache, report::NamedTuple)
    if length(X) != length(y)
        error("Inputs X and ouputs y have different lengths, ",
              "length(X) = $(length(X)), length(y) = $(length(y))")
    end
    crop_type = X.crop_type
    aquacrop_name = X.aquacrop_name
    soil_type = X.soil_type
    planting_density = X.planting_density
    weather_df = X.weather_df
    phenology_raw_df = X.phenology_raw_df
    years = X.years
    target_output = y

    println("\n")
    println("fitting years = $years")
    @show nrow(phenology_raw_df)

    # fit_type, sowing_phase, harvest_phase
    fit_type, sowing_phase, harvest_phase = get_fit_type_and_phases(crop_type)
    println("fit_type = $fit_type, sowing_phase = $sowing_phase, harvest_phase = $harvest_phase")

    # check sowing_phase, harvest_phase
    phases = unique(phenology_raw_df[!, :phase_id])
    if sowing_phase ∉ phases
        error("sowing_phase = $sowing_phase not found in phases = $phases")
    end
    if harvest_phase ∉ phases
        error("harvest_phase = $harvest_phase not found in phases = $phases")
    end

    # create a kw tuple with the additional information that we wish to pass to AquaCrop
    kw = (; crop_dict = Dict{String,Any}("PlantingDens" => planting_density))

    # 1. Calibrate Phenology
    println("running CropGrowthTutorial.calibrate_phenology_parameters()")
    @time crop_dict, phenology_df = CropGrowthTutorial.calibrate_phenology_parameters(
        phenology_raw_df, aquacrop_name, weather_df, sowing_phase, harvest_phase; kw...)
    if isnothing(phenology_df)
        println("Calibrating phenology parameters failed, not enough phenology data")
        error("Calibrating phenology parameters failed, not enough phenology data")
    end
    filter!(row -> row.year in years, phenology_df)
    @show nrow(phenology_df)
    println("phenology_df =")
    display(phenology_df)
    println("phenology_raw_df =")
    display(phenology_raw_df)

    # 2. Calibrate the canopy growth parameters
    println("running CropGrowthTutorial.calibrate_canopy_root_parameters!")
    @time CropGrowthTutorial.calibrate_canopy_root_parameters!(crop_dict, aquacrop_name);

    # 3. Calibrate biomass-yield parameters
    println("running CropGrowthTutorial.calibrate_biomass_yield_parameters!")
    @time CropGrowthTutorial.calibrate_biomass_yield_parameters!(
        crop_dict, aquacrop_name, soil_type, weather_df, fit_type, phenology_df, target_output;
        rng_seed=model.rng_seed,
        maxiter=model.maxiter,
    )
    
    # fit_description
    desc = Dict(
        :PhenologyCanopy => "phenology canopy calibration",
        :Biomass => "phenology canopy biomass calibration",
        :Yield => "phenology canopy yield calibration",
    )
    fit_description = desc[fit_type]

    # fitresult = crop
    fitresult = crop_dict
    cache = nothing
    report = (; fit_description, phenology_df)  # can put extra information here
    return fitresult, cache, report
end

# returns yhat (scalar)
function predict_one(crop_dict::Dict, weather_df::DataFrame, aquacrop_name::String,
                     soil_type::String, sowing_date::Date, fit_type::Symbol)
    cropfield, all_ok = CropGrowthTutorial.run_simulation(
        soil_type, aquacrop_name, sowing_date, weather_df;
        crop_dict
    )
    if !all_ok.logi
        return Inf
    else
        if fit_type == :Yield
            return cropfield.dayout[end,"Y(fresh)"].val
        elseif fit_type == :Biomass
            return cropfield.dayout[end,"Biomass"].val
        else
            error("Unknown fit_type $fit_type")
        end 
    end
end

function MMI.predict(model::AquaCropFit, fitted::Dict, X::CropInputs)
    # fit_type, sowing_phase, harvest_phase
    fit_type, sowing_phase, harvest_phase = get_fit_type_and_phases(X.crop_type)

    # sowing_dates
    phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(X.crop_type, X.station_name)
    sowing_dates = Dict{Int, Date}()  # year => sowing_date
    phenology_df = CropGrowthTutorial.process_crop_phenology(phenology_raw_df, sowing_phase, harvest_phase)
    for year in X.years
        df = filter(r -> r.year == year, phenology_df)
        if nrow(df) != 1
            @show phenology_df
            @show X.years
            error("Expected one row in dataframe, got $(nrow(df)) for year = $year")
        end
        sowing_dates[year] = df[1, :sowingdate]
    end

    return [predict_one(fitted, X.weather_df, X.aquacrop_name, X.soil_type, sowing_dates[year], fit_type)
            for year in X.years]
end



###################
# cross-validation
###################

function crossvalidate_fit(
    station_name::AbstractString,
    crop_type::AbstractString;
    years::Union{AbstractVector{Int}, Nothing}=nothing,
    folds::AbstractVector,
    repeats::Int=1,
    measures=DEFAULT_MEASURES,
    rng_seed::Int=rand(Int),
    maxiter::Int=DEFAULT_MAXITER,
)
    println()
    println("make_fit inputs: station_name = $station_name, crop_type = $crop_type")

    # aquacrop_name, planting_density
    crop_data_df = CropGrowthTutorial.get_crop_parameters()
    filter!(row -> row.crop_type == crop_type, crop_data_df)
    if nrow(crop_data_df) != 1
        error("Expected one matching row in crop_data_df for crop_type = $crop_type,",
              " found $(nrow(crop_data_df)) rows")
    end
    aquacrop_name = first(crop_data_df).aquacrop_name
    planting_density = first(crop_data_df).plantingdens  # TODO: units are num_seeds/ha ?

    # soil_type
    stations_df = CropGrowthTutorial.get_phenology_stations()
    filter!(row -> row.station_name == station_name, stations_df)
    if nrow(stations_df) != 1
        error("Expected one matching row in stations_df for station_name = $station_name,",
              " found $(nrow(stations_df)) rows")
    end
    soil_type = first(stations_df).soil_type

    # weather_df
    weather_df = CropGrowthTutorial.get_climate_data(station_name)
    if isnothing(weather_df)
        error("No climate data for station $station_name")
    end

    # yield_df
    yield_df = CropGrowthTutorial.get_yield_data(station_name)
    if isnothing(yield_df)
        error("No yield data for station $station_name")
    end
    # filter crop_type
    filter!(row -> row.crop_type == crop_type, yield_df)
    println("yield_df =")
    display(yield_df)

    # phenology_raw_df: crop phenology raw data for a given station
    phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
    if isnothing(phenology_raw_df)
        error("No phenology data for station $station_name and crop_type $crop_type")
    end

    # years
    if isnothing(years)
        # check for which years we have yield and phenology data
        years_yield_df = sort([parse(Int, colname) for colname in names(yield_df) if colname ∉ ("crop_type", "unit")])
        years_phenology_raw_df = sort(unique(phenology_raw_df.referenzjahr))
        years = intersect(years_yield_df, years_phenology_raw_df)
    end
    sort!(years)
    let allyears = collect(minimum(years):maximum(years))
        # TODO: duplicated
        years_yield_df = sort([parse(Int, colname) for colname in names(yield_df) if colname ∉ ("crop_type", "unit")])
        years_phenology_raw_df = sort(unique(phenology_raw_df.referenzjahr))
        if years != allyears
            println("Some years are missing\n" *
                "years = $years\n" *
                "missing years are = $(setdiff(allyears, years))\n" *
                "years_yield_df = $years_yield_df\n" *
                "years_phenology_raw_df = $years_phenology_raw_df\n")
        end
    end
    years::AbstractVector{Int}
    println("years = $years")

    # filter years from yield_df, which are stored as columns
    cols = ["crop_type", "unit", string.(years)...]
    yield_df = yield_df[:, cols]
    if nrow(yield_df) != 1
        error("Expected exactly one row in yield_df, crop_type = $crop_type at station_name = $station_name\n",
              "yield_df:\n$yield_df")
    end

    # yields
    yields = first(yield_df)

    # filter years from phenology_raw_df
    # We add one extra year, so that the phenology data for the
    # harvest date (which is in the following year) is included
    phenology_raw_df = filter(row -> row.referenzjahr in [years..., years[end]+1], phenology_raw_df)

    model = AquaCropFit(; rng_seed, maxiter)
    X = CropInputs(; station_name, crop_type, aquacrop_name, soil_type, planting_density,
                   weather_df, phenology_raw_df, years)
    # the yield/biomass data is in units of dt = 10 tons, so we divide by 10 to get tons
    # (actually it's dt/ha, and dividing by 10 gives ton/ha)
    y = Float64[getindex(yields, Symbol(year)) / 10 for year in years]

    println("soil_type = $soil_type, aquacrop_name = $aquacrop_name, planting_density = $planting_density")
    println("y = $y")
    println("folds = $folds")

    # Note: no data shuffling as data is time-ordered
    mach = MLJ.machine(model, X, y)
    #cv = MLJ.CV(; nfolds, shuffle=false)
    #cv = TimeSeriesCV(; nfolds)
    r = MLJ.evaluate!(
        mach;
        resampling=repeat(folds, repeats),  # we have to manually repeat the folds
        force = true,                       # needed because we have to manually repeat the folds
        measures,
        verbosity=1
    )
    return r
end

# Use CPU threads in MLJ
MLJ.default_resource(MLJ.CPUThreads())

# Available yield data for stations and crop types
#
# crop_type           Eichsfeld_years    Jena_years    Thüringer_Becken_years
# ---------           ---------------    ----------    ----------------------
# Silage maize        2009--2022         2000--2014    2014--2022
# Winter wheat        2009--2022         2000--2014    2014--2022
# Winter barley       2009--2022         2000--2014    2015--2022
# Winter rapeseed     2009--2022         2000--2014    2014--2022

function validate_all(; repeats::Int=1, measures=DEFAULT_MEASURES, rng_seed::Int=42, maxiter::Int=DEFAULT_MAXITER)
    rng = StableRNG(hash(rng_seed, UInt(0x1a2b3c4d5e6f7788)))  # salted hash

    inputs = Dict(
        "Eichsfeld" => Dict(
            # crop_type        years folds
            "silage_maize" => ([2009:2018..., 2020:2021...], [(1:10, 11:12)]),  # 2018: phenology data missing
            "winter_wheat" => ([2009:2018..., 2020:2021...], [(1:10, 11:12)]),  # 2018: phenology data missing
            "winter_barley" => ([], []),
            "winter_rape" => ([2009:2021...], [(1:10, 11:13)]),
        ),
        "Jena" => Dict(
            "silage_maize" => ([2000:2014...], [(1:10, 11:15)]),
            "winter_wheat" => ([2000:2013...], [(1:10, 11:14)]),  # 2014: problems with phenology data
            "winter_barley" => ([], []),  # sowing_phase not found
            "winter_rape" => ([2000:2013...], [(1:10, 11:14)]),  # 2014: problems with phenology data
        ),
        "Thuringer_Becken" => Dict(
            "silage_maize" => ([2014:2022...], [(1:6, 7:9)]),
            "winter_wheat" => ([2014, 2016:2021...], [(1:5, 6:7)]),  # 2015, 2022: problems, missing sowingdate
            "winter_barley" => ([2016:2021...], [(1:4, 5:6)]),  # 2014, 2015, 2022: missing data
            "winter_rape" => ([2014:2021...], [(1:6, 7:8)]),  # 2022: problems
        ),
    )
    rng_seeds = Dict(
        (station_name, crop_name) => rand(rng, Int)
        for station_name in keys(inputs) for crop_name in keys(inputs[station_name])
    )
    @show rng_seeds

    results = Dict{Tuple{String,String},Any}()
    fails = Any[]
    change_from_ref = Dict{Tuple{String,String}, Float64}()

    lk = ReentrantLock()
    @sync for station_name in sort(collect(keys(inputs)))
        crop_years_folds = inputs[station_name]
        println()
        @show station_name
        for (crop_name, (years, folds)) in crop_years_folds
            if isempty(folds) || isempty(years)
                @info "Skipping $crop_name, folds or years are empty..."
                continue
            end
            Threads.@spawn begin
                reference_params_path =
                    joinpath(pkgdir(CropGrowthTutorial), "data", "sims", "$(crop_name)_$(station_name).TOML")
                reference_params = TOML.parsefile(reference_params_path)["crop"]
                println()
                @show crop_name
                @show years
                @show folds
                for (f1, f2) in folds
                    @show years[f1], years[f2]
                end
                r = try
                    crossvalidate_fit(
                        station_name, crop_name;
                        folds, years, repeats, measures, maxiter,
                        rng_seed = rng_seeds[(station_name, crop_name)],
                    )
                catch err
                    @warn("Failed crossvalidation for $station_name and $crop_name",
                          exception=(err, catch_backtrace()))
                    lock(lk) do
                        push!(fails, (station_name, crop_name, years, folds))
                    end
                    nothing
                end
                changes = if !isnothing(r)
                    mean(mean_abs_rel_change(reference_params, nt.fitresult) for nt in r.fitted_params_per_fold)
                else
                    nothing
                end

                lock(lk) do
                    results[(station_name, crop_name)] = r
                    change_from_ref[(station_name, crop_name)] = changes
                end
            end
        end
    end

    # changes_df DataFrame
    changes_df = DataFrame([
        (station=a, crop=b, mean_abs_rel_change=c)
        for ((a,b), c) in change_from_ref if !isnothing(change_from_ref)
    ])
    sort!(changes_df, [:station, :crop])

    # results_perf_df DataFrame
    results_perf_df = DataFrame([
        let
            names = Symbol.(getindex.(split.(string.(perf_evals.measure), "("), 1))
            replace!(names,
                :RootMeanSquaredError => :RMSE,
                :LPLoss => :MAE,
                :MeanAbsoluteProportionalError => :MAPE,
                :RootMeanSquaredProportionalError => :RMSPE,
            )
            vals = perf_evals.measurement
            metrics = NamedTuple( (name => val for (name, val) in zip(names, vals)) )
            (; station_name, crop_name, metrics...)
        end
        for ((station_name, crop_name), perf_evals) in results if !isnothing(perf_evals)
    ])
    sort!(results_perf_df, [:station_name, :crop_name])

    return (; results, results_perf_df, changes_df, fails)
end

function show_years(station_name::AbstractString, crop_type::AbstractString)
    # yield_df
    yield_df = CropGrowthTutorial.get_yield_data(station_name)
    if isnothing(yield_df)
        error("No yield data for station $station_name")
    end
    # filter crop_type
    filter!(row -> row.crop_type == crop_type, yield_df)

    # phenology_raw_df: crop phenology raw data for a given station
    phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
    if isnothing(phenology_raw_df)
        error("No phenology data for station $station_name and crop_type $crop_type")
    end

    @show years_yield_df = sort([parse(Int, colname) for colname in names(yield_df) if colname ∉ ("crop_type", "unit")])
    @show years_phenology_raw_df = sort(unique(phenology_raw_df.referenzjahr))
    @show years = intersect(years_yield_df, years_phenology_raw_df)

    return yield_df, phenology_raw_df
end

function mean_abs_rel_change(d1::Dict{Tkey, T1}, d2::Dict{Tkey, T2}) where {Tkey, T1, T2}
    common_keys = intersect(keys(d1), keys(d2))
    r = 0.0
    if length(common_keys) == 0
        return 0.0
    end
    for k in common_keys
        r += abs((d2[k] - d1[k]) / d1[k])
    end
    return r / length(common_keys)
end
