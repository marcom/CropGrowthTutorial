###############
# DESCRIPTION #
###############
"""
Generate one multi-panel goodness-of-fit figure for multiple crops and all available
local-parameter stations.

For each crop/station pair with a calibrated TOML file in `data/sims`, the script:
- computes observed phenology durations from raw phenology + climate data,
- computes simulated phenology durations using the TOML crop parameters,
- compares observed vs simulated crop output (yield for non-silage, biomass for silage).

The resulting figure has 4 panels per crop:
- emergence (actual vs simulated days)
- begin flowering (actual vs simulated days)
- harvest (actual vs simulated days)
- output (actual vs simulated, in t/ha)

Usage:
    julia --project=. -t 4 scripts/plot_aquacrop_allregions_gof.jl

Year alignment and fallback caveats:
- Output-year handling intentionally differs by data availability:
  - if station phenology exists, output GOF is evaluated on phenology rows only;
  - if station phenology is missing, output GOF uses all available yield years.
- Yield lookup in the phenology path uses harvest year (`year(harvestdate)` when present).
- Fallback output simulation (no station phenology) uses a fixed sowing-year heuristic:
  - silage crops: sowing year = yield year;
  - non-silage crops: sowing year = yield year - 1.
- `DOY` means "day of year". Fallback sowing dates are reconstructed from cached reference
  sowing DOYs (median by year when available, otherwise crop default DOY).
- The non-silage fallback heuristic is pragmatic and may be imperfect for rare same-year
  sowing/harvest anomalies.
"""

# SETUP #
using DrWatson
@quickactivate "CropGrowthTutorial"

using AquaCrop
using CropGrowthTutorial
using CairoMakie
using DataFrames
using Dates
using Statistics
using Unitful

#############
# CONSTANTS #
#############
const PHENO_METRICS = [
    (:emergence, "days to emergence", :emergence_actualdays, :emergence_simulateddays),
    (:beginflowering, "days to begin flowering", :beginflowering_actualdays, :beginflowering_simulateddays),
    (:harvest, "days to harvest", :harvest_actualdays, :harvest_simulateddays),
]

const OUTPUT_METRIC = :output

const TOML_FIELDS = [
    "GDDCGC", "GDDaysToGermination", "GDDaysToFlowering", "GDDaysToMaxRooting",
    "GDDLengthFlowering", "GDDCDC", "RootMax", "GDDaysToHarvest", "GDDaysToFullCanopy",
    "GDtranspLow", "PlantingDens", "GDDaysToSenescence", "WP", "HI"
]

const DEFAULT_REGIONS = [
    "Bodensee",
    "Eichsfeld",
    "Hohenlohe",
    "Jena",
    "Oberrhein",
    "Thuringer_Becken",
]

const CROP_LABELS = Dict(
    "silage_maize" => "Silage maize",
    "winter_wheat" => "Winter wheat",
    "winter_barley" => "Winter barley",
    "winter_rape" => "Winter rape",
)

# Normalize known station-name variants used across data sources.
const STATION_ALIASES = Dict(
    "Thueringer_Becken" => "Thuringer_Becken",
)

const STATION_LABELS = Dict(
    "Thuringer_Becken" => "Thüringer Becken",
)

const GERMANY_PARAM_STATIONS = Set(["Bodensee", "Hohenlohe"])

const GERMANY_PARAM_FILES = Dict(
    "silage_maize" => "maize.toml",
    "winter_wheat" => "winter_wheat.toml",
    "winter_barley" => "winter_barley.toml",
    "winter_rape" => "winter_rape.toml",
)

const DEFAULT_SOWING_DOY = Dict(
    "silage_maize" => 120,
    "winter_wheat" => 285,
    "winter_barley" => 280,
    "winter_rape" => 245,
)

const SOWING_DOY_CACHE = Dict{String,Tuple{Dict{Int,Int},Int}}()

const STATION_COLORS = [
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", # green
    "#CC79A7", # magenta
    "#332288", # indigo
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#000000", # black
]

const STATION_MARKERS = [
    :circle,
    :rect,
    :utriangle,
    :diamond,
    :cross,
    :xcross,
    :star5,
    :hexagon,
]

const A4_TICKLABEL_SIZE = 27
const A4_AXISLABEL_SIZE = 31
const A4_TITLE_SIZE = 28
const A4_LEGEND_SIZE = 30

#############
# HELPERS   #
#############
function _keep_pairs(x, y)
    xx = Float64[]
    yy = Float64[]
    for (xi, yi) in zip(x, y)
        if ismissing(xi) || ismissing(yi)
            continue
        end
        xi_f = Float64(xi)
        yi_f = Float64(yi)
        if isfinite(xi_f) && isfinite(yi_f)
            push!(xx, xi_f)
            push!(yy, yi_f)
        end
    end
    return xx, yy
end

_crop_label(crop_type::AbstractString) = get(CROP_LABELS, crop_type, replace(crop_type, "_" => " "))
_is_silage(crop_type::AbstractString) = startswith(crop_type, "silage")
_output_metric_label(crop_type::AbstractString) = _is_silage(crop_type) ? "biomass [t/ha]" : "yield [t/ha]"
_ordered_unique(v::AbstractVector) = unique(v)

_canonical_station_name(station_name::AbstractString) = get(STATION_ALIASES, station_name, station_name)

function _station_name_candidates(station_name::AbstractString)
    canonical = _canonical_station_name(station_name)
    candidates = String[canonical]
    station_name != canonical && push!(candidates, station_name)
    for (alias, target) in STATION_ALIASES
        target == canonical && push!(candidates, alias)
    end
    return unique(candidates)
end

function _resolve_station_name(station_name::AbstractString, available_stations::AbstractSet{<:AbstractString})
    for candidate in _station_name_candidates(station_name)
        candidate in available_stations && return candidate
    end
    return nothing
end

_station_label(station_name::AbstractString) = begin
    canonical = _canonical_station_name(station_name)
    get(STATION_LABELS, canonical, replace(canonical, "_" => " "))
end

function _station_style_maps(station_set::AbstractVector{<:AbstractString})
    station_to_color = Dict{String,Any}()
    station_to_marker = Dict{String,Symbol}()
    for (i, st) in enumerate(station_set)
        station_to_color[st] = STATION_COLORS[mod1(i, length(STATION_COLORS))]
        station_to_marker[st] = STATION_MARKERS[mod1(i, length(STATION_MARKERS))]
    end
    return station_to_color, station_to_marker
end

function _load_crop_dict_from_toml(crop_type::AbstractString, station_name::AbstractString)
    for station_candidate in _station_name_candidates(station_name)
        toml_path = datadir("sims", crop_type * "_" * station_candidate * ".TOML")
        if !isfile(toml_path)
            continue
        end
        crop = AquaCrop.RepCrop()
        AquaCrop.load_gvars_from_toml!(crop, toml_path)
        return Dict{String,Any}(string(k) => getfield(crop, k) for k in fieldnames(typeof(crop)) if string(k) in TOML_FIELDS)
    end
    return nothing
end

function _has_local_toml(crop_type::AbstractString, station_name::AbstractString)
    for station_candidate in _station_name_candidates(station_name)
        toml_path = datadir("sims", crop_type * "_" * station_candidate * ".TOML")
        isfile(toml_path) && return true
    end
    return false
end

function _load_crop_dict_from_path(toml_path::AbstractString)
    isfile(toml_path) || return nothing
    crop = AquaCrop.RepCrop()
    AquaCrop.load_gvars_from_toml!(crop, toml_path)
    return Dict{String,Any}(string(k) => getfield(crop, k) for k in fieldnames(typeof(crop)) if string(k) in TOML_FIELDS)
end

function _load_crop_dict_from_germany_params(crop_type::AbstractString)
    param_file = get(GERMANY_PARAM_FILES, crop_type, nothing)
    isnothing(param_file) && return nothing
    return _load_crop_dict_from_path(datadir("germany-params", param_file))
end

function _load_crop_dict_for_station(crop_type::AbstractString, station_name::AbstractString)
    local_dict = _load_crop_dict_from_toml(crop_type, station_name)
    if local_dict !== nothing
        return local_dict, :local
    end
    canonical = _canonical_station_name(station_name)
    if canonical in GERMANY_PARAM_STATIONS
        germany_dict = _load_crop_dict_from_germany_params(crop_type)
        if germany_dict !== nothing
            return germany_dict, :germany
        end
    end
    return nothing, :missing
end

function _reference_sowing_doys(crop_type::AbstractString)
    if haskey(SOWING_DOY_CACHE, crop_type)
        return SOWING_DOY_CACHE[crop_type]
    end

    default_doy = get(DEFAULT_SOWING_DOY, crop_type, 120)
    raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type)
    raw_df === nothing && return (Dict{Int,Int}(), default_doy)

    harvest_phase = _is_silage(crop_type) ? 39 : 24
    phases = unique(raw_df[!, :phase_id])
    if !(10 in phases) || !(harvest_phase in phases)
        return (Dict{Int,Int}(), default_doy)
    end

    pheno_df = CropGrowthTutorial.process_crop_phenology(raw_df, 10, harvest_phase)
    if nrow(pheno_df) == 0
        return (Dict{Int,Int}(), default_doy)
    end

    per_year = Dict{Int,Vector{Int}}()
    all_doys = Int[]
    for row in eachrow(pheno_df)
        # Cache reference sowing timing as DOY, keyed by the sowing-season year.
        doy = dayofyear(row.sowingdate)
        push!(get!(per_year, row.year, Int[]), doy)
        push!(all_doys, doy)
    end
    if isempty(all_doys)
        return (Dict{Int,Int}(), default_doy)
    end

    year_to_doy = Dict{Int,Int}(year => round(Int, median(doys)) for (year, doys) in per_year)
    overall_doy = round(Int, median(all_doys))
    SOWING_DOY_CACHE[crop_type] = (year_to_doy, overall_doy)
    return SOWING_DOY_CACHE[crop_type]
end

function _fallback_sowingdate(crop_type::AbstractString, year::Int)
    # Rebuild a sowing date from year + reference DOY. The year passed here is decided by
    # _sowing_year_for_yield in the yield-only fallback path.
    year_to_doy, default_doy = _reference_sowing_doys(crop_type)
    doy = get(year_to_doy, year, default_doy)
    doy = clamp(doy, 1, dayofyear(Date(year, 12, 31)))
    return Date(year, 1, 1) + Day(doy - 1)
end

# In the phenology path, we align observed yield to harvest year when that date is available.
_yield_lookup_year(row) = !ismissing(row.harvestdate) ? year(row.harvestdate) : row.year
# Yield-only fallback heuristic: non-silage assumes sowing in previous calendar year.
_sowing_year_for_yield(crop_type::AbstractString, yield_year::Int) = _is_silage(crop_type) ? yield_year : (yield_year - 1)

function _phenology_eval_df_from_toml(
    crop_type::AbstractString,
    crop_name::AbstractString,
    station_name::AbstractString,
    hk_clim_df::DataFrame,
    crop_dict::Dict{String,Any},
)
    phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
    phenology_raw_df === nothing && return nothing

    sowing_phase = 10
    harvest_phase = _is_silage(crop_type) ? 39 : 24
    phases = unique(phenology_raw_df[!, :phase_id])
    if !(sowing_phase in phases) || !(harvest_phase in phases)
        return nothing
    end

    processed_pheno_df = CropGrowthTutorial.process_crop_phenology(phenology_raw_df, sowing_phase, harvest_phase)
    nrow(processed_pheno_df) > 0 || return nothing

    pheno_actual_df = CropGrowthTutorial.process_crop_phenology_actual_gdd(
        crop_name,
        processed_pheno_df,
        hk_clim_df;
        crop_dict,
    )
    nrow(pheno_actual_df) > 0 || return nothing

    pheno_simulated_df = CropGrowthTutorial.process_crop_phenology_simulated_gdd(
        crop_name,
        processed_pheno_df,
        hk_clim_df;
        crop_dict,
    )
    nrow(pheno_simulated_df) > 0 || return nothing

    return leftjoin(pheno_actual_df, pheno_simulated_df; on = :sowingdate)
end

function _yield_by_year(station_name::AbstractString, crop_type::AbstractString)
    hk_yield_df = CropGrowthTutorial.get_yield_data(station_name)
    hk_yield_df === nothing && return Dict{Int,Float64}()

    crop_yield_df = filter(row -> row.crop_type == crop_type, hk_yield_df)
    nrow(crop_yield_df) > 0 || return Dict{Int,Float64}()
    row = first(crop_yield_df)

    values = Dict{Int,Float64}()
    for col in names(row)
        year = tryparse(Int, String(col))
        isnothing(year) && continue
        raw_value = row[col]
        ismissing(raw_value) && continue
        values[year] = Float64(raw_value) / 10
    end
    return values
end

function _simulated_output(cropfield, crop_type::AbstractString)
    if _is_silage(crop_type)
        return ustrip(cropfield.dayout[end, "Biomass"])
    end
    return ustrip(cropfield.dayout[end, "Y(fresh)"])
end

function _phenology_and_yield_points(
    crop_type::AbstractString,
    crop_name::AbstractString,
    soil_type::AbstractString,
    station_name::AbstractString,
    crop_dict::Dict{String,Any},
)
    hk_clim_df = CropGrowthTutorial.get_climate_data(station_name)
    hk_clim_df === nothing && error("No climate data for station=$station_name")

    phenology_df = _phenology_eval_df_from_toml(crop_type, crop_name, station_name, hk_clim_df, crop_dict)
    yield_by_year = _yield_by_year(station_name, crop_type)

    metrics = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    if phenology_df !== nothing
        for (metric_sym, _, actual_col, simulated_col) in PHENO_METRICS
            x, y = _keep_pairs(phenology_df[!, actual_col], phenology_df[!, simulated_col])
            isempty(x) && continue
            metrics[metric_sym] = (x, y)
        end
    end

    if !isempty(yield_by_year)
        actual_output = Float64[]
        simulated_output = Union{Missing,Float64}[]

        if phenology_df === nothing
            # No station phenology: evaluate output GOF against all observed yield years using
            # fallback sowing dates inferred from year heuristic + reference/default DOY.
            for yield_year in sort(collect(keys(yield_by_year)))
                push!(actual_output, yield_by_year[yield_year])
                sowing_year = _sowing_year_for_yield(crop_type, yield_year)
                sowingdate = _fallback_sowingdate(crop_type, sowing_year)

                cropfield, all_ok = CropGrowthTutorial.run_simulation(
                    soil_type,
                    crop_name,
                    sowingdate,
                    hk_clim_df;
                    crop_dict,
                )
                if !all_ok.logi
                    push!(simulated_output, missing)
                    continue
                end
                push!(simulated_output, _simulated_output(cropfield, crop_type))
            end
        else
            # Station phenology exists: output GOF is restricted to phenology rows and their
            # corresponding yield years, which can exclude yield years not represented in
            # phenology observations.
            for row in eachrow(phenology_df)
                actual_val = get(yield_by_year, _yield_lookup_year(row), nothing)
                isnothing(actual_val) && continue
                push!(actual_output, actual_val)

                cropfield, all_ok = CropGrowthTutorial.run_simulation(
                    soil_type,
                    crop_name,
                    row.sowingdate,
                    hk_clim_df;
                    crop_dict,
                )
                if !all_ok.logi
                    push!(simulated_output, missing)
                    continue
                end
                push!(simulated_output, _simulated_output(cropfield, crop_type))
            end
        end

        x, y = _keep_pairs(actual_output, simulated_output)
        isempty(x) || (metrics[OUTPUT_METRIC] = (x, y))
    end

    return metrics
end

function _collect_all_data(crops::AbstractVector{<:AbstractString}, stations::AbstractVector{<:AbstractString})
    crop_data = CropGrowthTutorial.get_crop_parameters()
    crop_rows = Dict(row.crop_type => row for row in eachrow(crop_data))

    pheno_stations = CropGrowthTutorial.get_phenology_stations()
    available_station_names = Set(String(row.station_name) for row in eachrow(pheno_stations))
    station_soils = Dict(String(row.station_name) => (ismissing(row.soil_type) ? "missing" : String(row.soil_type)) for row in eachrow(pheno_stations))

    station_inputs = _ordered_unique(_canonical_station_name.(String.(stations)))
    all_data = Dict{String, Dict{String, Any}}()

    for crop_type in crops
        if !haskey(crop_rows, crop_type)
            @warn "Skipping $crop_type (not found in crop_data_general)"
            continue
        end

        crop_row = crop_rows[crop_type]
        crop_name = String(crop_row.aquacrop_name)
        crop_dict_base = Dict("PlantingDens" => crop_row.plantingdens)
        all_data[crop_type] = Dict{String, Any}()

        for station_input in station_inputs
            station_name = _resolve_station_name(station_input, available_station_names)
            if isnothing(station_name)
                @warn "Skipping station=$station_input (not in phenotype station list)"
                continue
            end

            soil_type = station_soils[station_name]
            toml_dict, param_source = _load_crop_dict_for_station(crop_type, station_input)
            if toml_dict === nothing
                @warn "No calibration parameters found for $(crop_type)/$(station_input)"
                continue
            end
            if param_source == :germany
                @info "Using germany-wide parameters for $(crop_type)/$(station_input)"
            end
            if !haskey(toml_dict, "PlantingDens")
                toml_dict["PlantingDens"] = crop_dict_base["PlantingDens"]
            end
            if !_is_silage(crop_type) && !haskey(toml_dict, "HI")
                @warn "Missing HI override for $(crop_type)/$(station_input); simulation will use crop default HI"
            end

            try
                station_metrics = _phenology_and_yield_points(
                    crop_type,
                    crop_name,
                    soil_type,
                    station_name,
                    toml_dict,
                )
                if isempty(station_metrics)
                    @warn "No comparable observations for $(crop_type)/$(station_name)"
                    continue
                end
                all_data[crop_type][_canonical_station_name(station_name)] = station_metrics
            catch err
                @warn "Skipping $(crop_type)/$(station_name): $err"
            end
        end
    end

    return all_data
end

function _stations_with_toml(crops::AbstractVector{<:AbstractString})
    toml_files = readdir(datadir("sims"))
    stations = String[]
    crop_set = Set(crops)

    for f in toml_files
        endswith(lowercase(f), ".toml") || continue
        stem = splitext(f)[1]

        for crop in crop_set
            prefix = string(crop, "_")
            if startswith(stem, prefix)
                station = stem[length(prefix)+1:end]
                if !isempty(station)
                    push!(stations, _canonical_station_name(station))
                end
                break
            end
        end
    end

    return sort(unique(stations))
end

function _validate_required_parameter_files(crops::AbstractVector{<:AbstractString}, regions::AbstractVector{<:AbstractString})
    missing_entries = String[]
    missing_files = Set{String}()

    for crop_type in crops
        germany_file = get(GERMANY_PARAM_FILES, crop_type, nothing)
        for region in regions
            station_name = _canonical_station_name(String(region))
            _has_local_toml(crop_type, station_name) && continue
            !(station_name in GERMANY_PARAM_STATIONS) && continue

            if isnothing(germany_file)
                push!(missing_entries, "$(crop_type)/$(station_name) -> missing GERMANY_PARAM_FILES entry")
                continue
            end

            required_file = datadir("germany-params", germany_file)
            if !isfile(required_file)
                push!(missing_entries, "$(crop_type)/$(station_name) -> $(required_file)")
                push!(missing_files, required_file)
            end
        end
    end

    isempty(missing_entries) && return nothing

    missing_file_lines = join(sort(collect(missing_files)), "\n  - ")
    missing_pair_lines = join(sort(unique(missing_entries)), "\n  - ")
    msg = """
Missing required germany fallback parameter files.
Selected regions include stations that require germany-wide parameters when no local TOML exists.

Missing files:
  - $(isempty(missing_file_lines) ? "n/a (mapping missing instead)" : missing_file_lines)

Affected crop/station pairs:
  - $(missing_pair_lines)
"""
    error(msg)
end

################
# PLOTTING     #
################
function _panel_specs(crop_type::AbstractString)
    return [
        (:emergence, "days to emergence"),
        (:beginflowering, "days to begin flowering"),
        (:harvest, "days to harvest"),
        (OUTPUT_METRIC, _output_metric_label(crop_type)),
    ]
end

function make_figure(
    all_data;
    crops::AbstractVector{<:AbstractString},
    station_order::AbstractVector{<:AbstractString},
    out_svg::AbstractString,
    out_png::AbstractString,
)
    station_lists = [collect(keys(v)) for v in values(all_data) if !isempty(v)]
    isempty(station_lists) && error("No data found for selected crops/stations")
    station_set = _ordered_unique(_canonical_station_name.(String.(station_order)))
    isempty(station_set) && error("No stations selected for plotting")

    station_to_color, station_to_marker = _station_style_maps(station_set)

    f = Figure(size = (2200, 1700))
    for (ri, crop_type) in enumerate(crops)
        for (ci, (metric_sym, metric_label)) in enumerate(_panel_specs(crop_type))
            ax = Axis(
                f[ri, ci],
                title = string(_crop_label(crop_type), " — ", metric_label),
                xlabel = "actual",
                ylabel = "simulated",
                aspect = DataAspect(),
                titlesize = A4_TITLE_SIZE,
                xlabelsize = A4_AXISLABEL_SIZE,
                ylabelsize = A4_AXISLABEL_SIZE,
                xticklabelsize = A4_TICKLABEL_SIZE,
                yticklabelsize = A4_TICKLABEL_SIZE,
            )

            panel_x = Float64[]
            panel_y = Float64[]
            crop_metrics = get(all_data, crop_type, Dict{String, Any}())

            for station_name in station_set
                st_data = get(crop_metrics, station_name, nothing)
                st_data === nothing && continue
                haskey(st_data, metric_sym) || continue

                x, y = st_data[metric_sym]
                isempty(x) && continue

                append!(panel_x, x)
                append!(panel_y, y)
                scatter!(
                    ax,
                    x,
                    y;
                    color = station_to_color[station_name],
                    marker = station_to_marker[station_name],
                    markersize = 9,
                    strokecolor = :black,
                    strokewidth = 0.4,
                    label = station_name,
                )
            end

            if !isempty(panel_x)
                lo = min(minimum(panel_x), minimum(panel_y))
                hi = max(maximum(panel_x), maximum(panel_y))
                if lo == hi
                    lo -= 1
                    hi += 1
                end
                lines!(ax, [lo, hi], [lo, hi]; color = :tomato, linestyle = :dash, linewidth = 1.2)
                xlims!(ax, lo, hi)
                ylims!(ax, lo, hi)
            end
        end
    end

    legend_elems = [
        MarkerElement(
            color = station_to_color[s],
            marker = station_to_marker[s],
            markersize = 10,
            strokecolor = :black,
            strokewidth = 0.4,
        ) for s in station_set
    ]
    legend_labels = [_station_label(s) for s in station_set]

    last_row = length(crops) + 1
    ncols = length(_panel_specs(first(crops)))
    f[last_row, 1:ncols] = Legend(
        f,
        legend_elems,
        legend_labels;
        labelsize = A4_LEGEND_SIZE,
        orientation = :horizontal,
        tellwidth = false,
        tellheight = true,
        framevisible = false,
    )

    save(out_svg, f)
    save(out_png, f)
    return f
end

############
# RUNTIME  #
############
function main(; crops = ["silage_maize", "winter_wheat", "winter_barley", "winter_rape"], regions = nothing)
    regions_ = isnothing(regions) ? copy(DEFAULT_REGIONS) : _ordered_unique(_canonical_station_name.(String.(regions)))
    _validate_required_parameter_files(crops, regions_)
    data = _collect_all_data(crops, regions_)

    out_svg = plotsdir("aquacrop_reliability_all_regions.svg")
    out_png = plotsdir("aquacrop_reliability_all_regions.png")
    mkpath(dirname(out_svg))

    return make_figure(
        data;
        crops = crops,
        station_order = regions_,
        out_svg = out_svg,
        out_png = out_png,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
