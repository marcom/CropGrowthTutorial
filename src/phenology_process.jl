#####################################################
# functions to process and calibrate phenology data #
#####################################################
"""
    df = process_crop_phenology(raw_phenology_df::AbsDataFrame, sowing_phase::Int, harvest_phase::Int)

returns a DataFrame with the dates of the `sowing_phase` and `harvest_phase` and other important phenology phases 
"""
function process_crop_phenology(raw_phenology_df::AbstractDataFrame, sowing_phase::Int, harvest_phase::Int)
    # this phases are from data/exp_raw/phase_descriptions.csv
    additional_phases = Dict(
        :emergencedate=>[12, 16], 
        :beginfloweringdate=>[5, 18, 17], 
        :endfloweringdate=>[7, 6, 41]) 
        # :senescensedate=>[21, 31, 32, 63, 64, 22])

    dateformat = DateFormat("yyyymmdd")
    make_date = (x) -> Date(string(x), dateformat)
    df = copy(raw_phenology_df)
    df[!,:date] = make_date.(df.eintrittsdatum)

    sort!(df, :date)
    sow_i = 1
    harv_i = 1
    full_cols = Dict(:sowingdate => Date[],
                     :harvestdate => Date[],
                     :daystoharvest => Int[],
                     :year => Int[],
                     :stations_id => Int[])
    other_cols = Dict(key=>Union{Missing,Date}[] for key in keys(additional_phases))
    outdf = DataFrame(merge(full_cols, other_cols))
    logi = true
    gb = groupby(df, :phase_id)
    phases_data = Dict(key => Dict(additional_phases[key][ii] =>  Vector(get(gb, (additional_phases[key][ii],), DataFrame(date=[]))[!,:date]) for ii in eachindex( additional_phases[key])) for key in keys(additional_phases))


    while logi
        if gb[(harvest_phase,)][harv_i,:date] <= gb[(sowing_phase,)][sow_i,:date]
            harv_i += 1
        elseif gb[(harvest_phase,)][harv_i,:date] >= gb[(sowing_phase,)][sow_i,:date] + Year(1)
            sow_i += 1
        else
            sowingdate = gb[(sowing_phase,)][sow_i,:date]
            harvestdate = gb[(harvest_phase,)][harv_i,:date]
            row = Dict(
                :sowingdate => sowingdate,
                :harvestdate => harvestdate,
                :daystoharvest => Day(harvestdate - sowingdate).value,
                :year => gb[(sowing_phase,)][sow_i,:referenzjahr],
                :stations_id => gb[(sowing_phase,)][sow_i,:stations_id],
            )
            for key in keys(additional_phases)
                for ii in eachindex(additional_phases[key])
                    if ismissing(get(row, key, missing))
                        i = findfirst(x->sowingdate<x<harvestdate, phases_data[key][additional_phases[key][ii]])
                        if isnothing(i)
                            phasedate = missing
                        else
                            phasedate = popat!(phases_data[key][additional_phases[key][ii]], i)
                        end
                        row[key] = phasedate
                    end
                end
            end


            push!(outdf, row)
            sow_i += 1
            harv_i += 1
        end
        if (sow_i>size(gb[(sowing_phase,)],1)) || (harv_i>size(gb[(harvest_phase,)],1))
            logi = false
        end
    end

   return outdf
end

"""
    df = process_crop_phenology_actual_gdd(crop_name, phenology_df, hk_clim_df; kw...)

returns a DataFrame with the duration in days and gddays of some important phenology phases using the actual data
"""
function process_crop_phenology_actual_gdd(crop_name, phenology_df, hk_clim_df; kw...)
    df = DataFrame(
                    sowingdate = Date[],
                    harvestdate = Date[],
                    year = Int[],

                    harvest_actualdays = Union{Missing,Int}[],
                    harvest_actualgdd = Union{Missing,Float64}[],
    
                    beginflowering_actualdays = Union{Missing,Int}[],
                    beginflowering_actualgdd = Union{Missing,Float64}[],

                    endflowering_actualdays = Union{Missing,Int}[],
                    endflowering_actualgdd = Union{Missing,Float64}[],

                    emergence_actualdays = Union{Missing,Int}[],
                    emergence_actualgdd = Union{Missing,Float64}[],
                    )

    # setup the crop since we need some parameters of it
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux = haskey(kw, :crop_dict) ? kw[:crop_dict] : nothing)
    
    for row in eachrow(phenology_df)
        gdd_df, all_ok = run_cropgdd(crop, row.sowingdate, row.harvestdate, hk_clim_df; kw...);
        if !all_ok.logi
            continue
        end
        res = find_actual_days_gdd(gdd_df, row)
        res[:year] = row.year
        push!(df, res)
    end

    return df
end

function find_actual_days_gdd(crop_gdd, row)
    res = Dict()

    res[:sowingdate] = row.sowingdate
    res[:harvestdate] = row.harvestdate

    for phase in [:harvestdate, :beginfloweringdate, :endfloweringdate, :emergencedate]
        gather_gdd_for_phase_actual!(res, row, phase, res[:sowingdate], crop_gdd)
    end

    return res
end

function gather_gdd_for_phase_actual!(res, row, phase, sowingdate, crop_gdd)
    sy_day = Symbol(replace(String(phase), "date" => "_actualdays"))
    sy_gdd = Symbol(replace(String(phase), "date" => "_actualgdd"))

    val = row[phase] - sowingdate
    if !ismissing(val)
        res[sy_day] = val.value
        res[sy_gdd] = find_day_gdd_for_date(row[phase], crop_gdd)
    else
        res[sy_day] = missing
        res[sy_gdd] = missing
    end

    return nothing 
end

function find_day_gdd_for_date(date, crop_gdd)
    rowi_data = filter(row->row.date==date, crop_gdd)
    return rowi_data[1, "GDD"]
end

"""
    df = process_crop_phenology_simulated_gdd(crop_name, phenology_df, hk_clim_df; kw...)

returns a DataFrame with the duration in days and gddays of some important phenology phases using the simulation
"""
function process_crop_phenology_simulated_gdd(crop_name, phenology_df, hk_clim_df; kw...)
    df = DataFrame(
                    sowingdate = Date[],

                    harvest_simulateddays = Union{Missing,Int}[],
                    harvest_simulatedgdd = Union{Missing,Float64}[],
    
                    beginflowering_simulateddays = Union{Missing,Int}[],
                    beginflowering_simulatedgdd = Union{Missing,Float64}[],

                    endflowering_simulateddays = Union{Missing,Int}[],
                    endflowering_simulatedgdd = Union{Missing,Float64}[],

                    emergence_simulateddays = Union{Missing,Int}[],
                    emergence_simulatedgdd = Union{Missing,Float64}[]
                    )
                    
    # setup the crop since we need some parameters of it
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux = haskey(kw, :crop_dict) ? kw[:crop_dict] : nothing)

    for row in eachrow(phenology_df)
        gdd_df, all_ok = run_cropgdd(crop, row.sowingdate, row.sowingdate + days_after_sowing, hk_clim_df; kw...);
        if !all_ok.logi
            continue
        end

        res = find_simulated_days_gdd(gdd_df, crop)
        res[:sowingdate] = row.sowingdate
        push!(df, res)
    end

    return df
end

function find_simulated_days_gdd(gdd_df, crop)
    res = Dict()

    sowingdate = gdd_df[1,:date]
    res[:sowingdate] = sowingdate 
    max_size = size(gdd_df,1)

    _, i = findmin( x -> abs(x-crop.GDDaysToHarvest), gdd_df[:,:GDD])
    if i<max_size
        res[:harvest_simulateddays] = (gdd_df[i, :date] - sowingdate).value
        res[:harvest_simulatedgdd] = gdd_df[i, :GDD]
    else
        res[:harvest_simulateddays] = missing
        res[:harvest_simulatedgdd] = missing
    end

    _, i = findmin( x -> abs(x-crop.GDDaysToFlowering), gdd_df[:,:GDD])
    if i<max_size
        res[:beginflowering_simulateddays] = (gdd_df[i, :date] - sowingdate).value
        res[:beginflowering_simulatedgdd] = gdd_df[i, :GDD]
    else
        res[:beginflowering_simulateddays] = missing
        res[:beginflowering_simulatedgdd] = missing
    end


    _, i = findmin( x -> abs(x-crop.GDDaysToFlowering-crop.GDDLengthFlowering), gdd_df[:,:GDD])
    if i<max_size
        res[:endflowering_simulateddays] = (gdd_df[i, :date] - sowingdate).value
        res[:endflowering_simulatedgdd] = gdd_df[i, :GDD]
    else
        res[:endflowering_simulateddays] = missing
        res[:endflowering_simulatedgdd] = missing
    end

    _, i = findmin( x -> abs(x-crop.GDDaysToGermination), gdd_df[:,:GDD])
    if i<max_size
        res[:emergence_simulateddays] = (gdd_df[i, :date] - sowingdate).value
        res[:emergence_simulatedgdd] = gdd_df[i, :GDD]
    else
        res[:emergence_simulateddays] = missing
        res[:emergence_simulatedgdd] = missing
    end

    return res
end

"""
    crop_dict, df = calibrate_phenology_parameters(raw_phenology_df, crop_name, hk_clim_df, sowing_phase, harvest_phase, method=:median; kw...)

returns a crop_dict with the calibrated crop parameters related to phenology stages
"""
function calibrate_phenology_parameters(raw_phenology_df, crop_name, hk_clim_df, sowing_phase, harvest_phase, method=:median; kw...)
    # check if we hace a crop_dict
    if haskey(kw, :crop_dict)
        crop_dict = kw[:crop_dict]
    else
        crop_dict = Dict()
    end

    pheno_df = process_crop_phenology(raw_phenology_df, sowing_phase, harvest_phase)
    if size(pheno_df,1)<minimal_pheno_data
        # not enough data for statistical fit, early return
        @warn "Not enough pheno data, size(pheno_df,1) == $(size(pheno_df,1)) < minimal_pheno_data == $minimal_pheno_data"
        display(pheno_df)
        return crop_dict, nothing
    end
    pheno_actual_df = process_crop_phenology_actual_gdd(crop_name, pheno_df, hk_clim_df; kw...)
     
    # setup the crop since we need some parameters of it
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux = haskey(kw, :crop_dict) ? kw[:crop_dict] : nothing)

    # use method to choose function to calibrate
    if method == :median
        ff = median
    elseif method == :mean
        ff = mean
    else
        ff = median
    end

    
    pars = ["GDDaysToHarvest", "GDDaysToFlowering", "GDDaysToGermination"]
    cols = [:harvest_actualgdd, :beginflowering_actualgdd, :emergence_actualgdd]
    for (par_, col) in zip(pars, cols)
        v = pheno_actual_df[!,col]
        crop_dict[par_] = _select_val( ff(v), getfield(crop, Symbol(par_)) )
    end

    v = pheno_actual_df[:,:endflowering_actualgdd]
    par_ = "GDDLengthFlowering"
    crop_dict[par_] = _select_val( ff(v) - crop_dict["GDDaysToFlowering"], getfield(crop, Symbol(par_)) ) 

    # create a kw tuple with the additional information that we wish to pass to AquaCrop
    # consider the median of the actual gdd distribution for each phenology phase
    kwargs = (
            crop_dict = crop_dict, 
         );

    # get the days and growing degree days for each phenology phase from the simulated data using the additional information
    pheno_simulated_df = process_crop_phenology_simulated_gdd(crop_name, pheno_df, hk_clim_df; kwargs...)

    # compare the day difference between the actual data and the simulated data
    df = leftjoin(pheno_actual_df, pheno_simulated_df; on=:sowingdate);

    return crop_dict, df
end

function _select_val(x, y)
    a = round(Int, coalesce(x, y))
    if a > 0 
        return a
    else
        return y
    end
end
