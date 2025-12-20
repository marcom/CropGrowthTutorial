###############
# DESCRIPTION #
###############
"""
    We make the crop calibration for all the crops in one of the stations

    Note:
    1. Some stations do not have enough phenology data to calibrate the phenology
    1. Some stations do not have yield/biomass data

    If we only have full data for a given crop in a given station we do the full calibration
    If we only have enough data for the phenology calibration we only do this
    If we do not have enough data for anything we skip the crop

    Results of crop calibration are shown has figures in  `plots/` file
"""
#########
# SETUP #
#########
using DrWatson
@quickactivate "CropGrowthTutorial"

# using Revise
using AquaCrop
using CropGrowthTutorial
using Dates
using DataFrames
using Unitful
using CairoMakie
using Plots
using Infiltrator

Plots.gr()
Plots.default(
    fontfamily="helvetica",
    guidefont=font(12),   # axes titles
    tickfont=font(10),   # tick labels
    legendfont=font(10),    # legend text
    guidefontfamily="helvetica",
    fontfamily_subplot="helvetica",
    legend_font_family="helvetica",
    legend_title_font_family="helvetica",
    plot_titlefontfamily="helvetica",
    titlefontfamily="helvetica",
    xguidefontfamily="helvetica",
    yguidefontfamily="helvetica",
    xtickfontfamilya="helvetica",
    ytickfontfamilya="helvetica",
)


########
# MAIN #
########
function main(station_=nothing)
    # DFAULT STATION #
    default_station = "Jena"
    plottype = CropGrowthTutorial.MakiePlotOption()

    ## READ DATA ##

    # Stations information
    stations_df = CropGrowthTutorial.get_phenology_stations()
    # Crop type general information
    general_crop_data_df = CropGrowthTutorial.get_crop_parameters()

    ddd = Dict()


    if isnothing(station_)
        println("\n\nWARNING: no station name given, selecting default station ", default_station)
        station_ = default_station
    elseif !(station_ in stations_df.station_name)
        println("\n\nWARNING: station name ", station_, " not found in satations_df, selecting default station ", default_station)
        station_ = default_station
    end


    # find given station
    row = filter(x -> x.station_name == station_, stations_df)[1, :]
    station_name = row.station_name
    soil_type = ismissing(row.soil_type) ? "missing" : row.soil_type
    fit_type = :unknown
    println("\n\nCCCCCCCCCCCCCCCCCCCCCCC\n", "station_name ", station_name)

    # get climate data for given station
    hk_clim_df = CropGrowthTutorial.get_climate_data(station_name)
    if isnothing(hk_clim_df)
        println("WARNING: No climate data for station ", station_name, " going to next station")
        return nothing
    end
    # get yield data for given station
    hk_yield_df = CropGrowthTutorial.get_yield_data(station_name)
    if isnothing(hk_yield_df)
        println("WARNING: No yield data for station ", station_name, " we will only make phenolgy and canopy fit")
        fit_type = :PhenologyCanopy
    end

    # for each crop
    for roww in eachrow(general_crop_data_df)
        crop_type = roww.crop_type
        crop_name = roww.aquacrop_name
        plantingdens = roww.plantingdens
        println("\n###\n", "crop_type ", crop_type)

        # Crop's phenology raw data for a given station
        phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
        if isnothing(phenology_raw_df)
            println("WARNING: No phenology data for station ", station_name,
                " in crop_type ", crop_type,
                " going to next crop")
            continue
        end


        ## FIT TYPE ##
        sowing_phase = 10
        if startswith(crop_type, "silage")
            if fit_type != :PhenologyCanopy
                println("INFO: crop_type is silage, we will make biomass fit")
                fit_type = :Biomass
            end
            harvest_phase = 39
        else
            if fit_type != :PhenologyCanopy
                println("INFO: crop_type is not silage, we will make yield fit")
                fit_type = :Yield
            end
            harvest_phase = 24
        end
        vv = unique(phenology_raw_df[!, :phase_id])
        if (!(sowing_phase in vv) || !(harvest_phase in vv))
            println("WARNING: No enough phenology data for station ", station_name,
                " in crop_type ", crop_type,
                " going to next crop")
            continue
        end


        # create a kw tuple with the additional information that we wish to pass to AquaCrop
        kw = (crop_dict=Dict{String,Any}("PlantingDens" => plantingdens),)


        ## CALIBRATION ##
        if fit_type == :unknown
            println("WARNING: Unknown fit type selected for station ", station_name,
                " in crop_type ", crop_type,
                " going to next crop")
            continue
        else
            # 1. Calibrate Phenology
            crop_dict, phenology_df = CropGrowthTutorial.calibrate_phenology_parameters(phenology_raw_df, crop_name, hk_clim_df,
                sowing_phase, harvest_phase; kw...)
            if isnothing(phenology_df)
                println("WARNING: No enough phenology data for station ", station_name,
                    " in crop_type ", crop_type,
                    " going to next crop")
                continue
            end

            ddd[crop_type] = Dict()

            # 2. Calibrate the canopy growth parameters
            CropGrowthTutorial.calibrate_canopy_root_parameters!(crop_dict, crop_name)

            if fit_type == :PhenologyCanopy
                fit_description = "phenology canopy calibration"
            else
                # gather the yield/biomass data
                # NOTE: we assume we have this data, otherwise we should check and if missing go back to :PhenologyCanopy fit
                yield_df = filter(row -> row.crop_type == crop_type, hk_yield_df)
                target_output = []
                _phenology_df = copy(phenology_df)
                inds = Int[]
                if size(yield_df, 1) == 0
                    println("WARNING: No yield data for crop_type=" * crop_type * "on station_name=" *
                            station_name * " going to next crop")
                    continue
                end

                for (_i, rowi) in enumerate(eachrow(phenology_df))
                    yeari = rowi.year
                    crop_date = rowi.sowingdate
                    if (hk_clim_df[1, "date"] < crop_date < hk_clim_df[end, "date"]) && (string(yeari) in names(yield_df)) && !(ismissing(yield_df[1, string(yeari)]))
                        append!(target_output, yield_df[1, string(yeari)] / 10)
                    else
                        push!(inds, _i)
                    end
                end
                if length(inds) > 0
                    deleteat!(_phenology_df, inds)
                end
                if length(target_output) == 0
                    println("WARNING: Not enough yield data for crop_type=" * crop_type * "on station_name=" *
                            station_name * " going to next crop")
                    continue
                end

                # 3. Calibrate biomass-yield parameters
                CropGrowthTutorial.calibrate_biomass_yield_parameters!(crop_dict, crop_name, soil_type, hk_clim_df, fit_type, _phenology_df, target_output)

                if fit_type == :Biomass
                    fit_description = "phenology canopy biomass calibration"
                elseif fit_type == :Yield
                    fit_description = "phenology canopy yield calibration"
                end
            end

            ## PLOT RESULTS ##
            # plot of days to harvest actual and simulated
            f = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.harvest_actualdays, _phenology_df.harvest_simulateddays, crop_type, station_name, "days to harvest")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages", crop_type * "_" * station_name * "harvest_correlation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "harvest_correlation.svg"), f)
                Makie.save(plotsdir(crop_type * "_" * station_name * "harvest_correlation.png"), f)
            end
            ddd[crop_type]["harvest_actual"] = _phenology_df.harvest_actualdays
            ddd[crop_type]["harvest_simulated"] = _phenology_df.harvest_simulateddays

            f = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.emergence_actualdays, _phenology_df.emergence_simulateddays, crop_type, station_name, "days to emergence")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages", crop_type * "_" * station_name * "emergence_correlation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "emergence_correlation.svg"), f)
                Makie.save(plotsdir(crop_type * "_" * station_name * "emergence_correlation.png"), f)
            end
            ddd[crop_type]["emergence_actual"] = _phenology_df.emergence_actualdays
            ddd[crop_type]["emergence_simulated"] = _phenology_df.emergence_simulateddays

            f = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.beginflowering_actualdays, _phenology_df.beginflowering_simulateddays, crop_type, station_name, "days to begin flowering")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages", crop_type * "_" * station_name * "beginflowering_correlation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "beginflowering_correlation.svg"), f)
                Makie.save(plotsdir(crop_type * "_" * station_name * "beginflowering_correlation.png"), f)
            end
            ddd[crop_type]["beginflowering_actual"] = _phenology_df.beginflowering_actualdays
            ddd[crop_type]["beginflowering_simulated"] = _phenology_df.beginflowering_simulateddays

            # f = CropGrowthTutorial.plot_GDD_stats_years(_phenology_df, station_name)
            # Makie.save(plotsdir(crop_type*"_"*station_name*"_phenology.png"), f)

            # plot of crop growth simulation
            i = CropGrowthTutorial.find_closest_to_median(target_output)
            sowingdate = _phenology_df[i, "sowingdate"]
            endday = _phenology_df[i, "harvestdate"]
            kw = (
                crop_dict=crop_dict,
                end_day=endday
            )
            if !ismissing(_phenology_df[i, "emergence_actualdays"])
                emergence_day = sowingdate + Day(_phenology_df[i, "emergence_actualdays"])
                kw = merge(kw, (emergence_day=emergence_day,))
            end
            if !ismissing(_phenology_df[i, "beginflowering_actualdays"])
                beginflowering_day = sowingdate + Day(_phenology_df[i, "beginflowering_actualdays"])
                kw = merge(kw, (beginflowering_day=beginflowering_day,))
            end
            if !ismissing(_phenology_df[i, "endflowering_actualdays"])
                endflowering_day = sowingdate + Day(_phenology_df[i, "endflowering_actualdays"])
                kw = merge(kw, (endflowering_day=endflowering_day,))
            end


            cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df; kw...)
            if all_ok.logi == false
                println("\nWARNING: bad simulation result on ", sowingdate)
                println("  ", all_ok.msg)
                continue
            end
            f = CropGrowthTutorial.plot_daily_stuff_one_year(plottype, cropfield, crop_type, soil_type; kw...)
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages", crop_type * "_" * station_name * "_cropevolution.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "_cropevolution.svg"), f)
                Makie.save(plotsdir(crop_type * "_" * station_name * "_cropevolution.png"), f)
            end

            # plot yield results
            if fit_type in [:Biomass, :Yield]
                years = _phenology_df[:, "year"]
                simulated_yield = []
                simulated_yield_baseline = []
                for sowingdate in _phenology_df[:, "sowingdate"]
                    # fited yield
                    cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df; kw...)
                    if all_ok.logi == false
                        println("\nWARNING: bad simulation result on ", sowingdate, " for fited yield")
                        println("  ", all_ok.msg)
                        append!(simulated_yield, missing)
                        continue
                    end
                    if fit_type == :Biomass
                        append!(simulated_yield, ustrip(cropfield.dayout[end, "Biomass"]))
                    elseif fit_type == :Yield
                        append!(simulated_yield, ustrip(cropfield.dayout[end, "Y(fresh)"]))
                    end

                    # base line yield
                    cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df)
                    if all_ok.logi == false
                        println("\nWARNING: bad simulation result on ", sowingdate, " for baseline yield")
                        println("  ", all_ok.msg)
                        append!(simulated_yield_baseline, missing)
                        continue
                    end
                    if fit_type == :Biomass
                        append!(simulated_yield_baseline, ustrip(cropfield.dayout[end, "Biomass"]))
                    elseif fit_type == :Yield
                        append!(simulated_yield_baseline, ustrip(cropfield.dayout[end, "Y(fresh)"]))
                    end
                end

                dats = [target_output, simulated_yield, simulated_yield_baseline]
                cols = ["actual yield", "fitted yield", "baseline yield"]

                f = CropGrowthTutorial.plot_correlation(plottype, target_output, simulated_yield, crop_type, station_name, "yield data")
                if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                    Plots.savefig(f, plotsdir("svgimages", crop_type * "_" * station_name * "yield_correlation.svg"))
                else
                    Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "yield_correlation.svg"), f)
                    Makie.save(plotsdir(crop_type * "_" * station_name * "yield_correlation.png"), f)
                end
                ddd[crop_type]["yield_actual"] = target_output
                ddd[crop_type]["yield_simulated"] = simulated_yield

                f = CropGrowthTutorial.plot_yearly_data(dats, cols, years, crop_type, soil_type, station_name)
                Makie.save(plotsdir(crop_type * "_" * station_name * "_yearlydata.png"), f)
            end
            println("INFO: saved plots of crop_type " * crop_type * " from station_name " * station_name)
        end
    end
    return ddd
end

function write_fited_crop(crop_type, soil_type, crop_name, crop_dict, station_name, fit_description)
    # file where to save result
    file_name = crop_type * "_" * station_name * ".TOML"
    full_filename = datadir("sims", file_name)

    # make crop
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux=crop_dict)

    # make header
    head = "Crop fited: crop_type=" * crop_type * " station_name=" * station_name *
           " soil_type=" * soil_type * " " * fit_description

    # save crop
    AquaCrop.save_crop(full_filename, crop, head)
    return nothing
end


function plotting(ddd, stations, crop_type)
    # 4 panels in a 2×2
    stuffs = ["emergence", "beginflowering", "harvest", "yield"]

    # Map each stuff → (actual, simulated) keys
    select_stuff(stuff) = begin
        if stuff == "yield"
            "yield_actual", "yield_simulated"
        elseif stuff == "beginflowering"
            "beginflowering_actual", "beginflowering_simulated"
        elseif stuff == "harvest"
            "harvest_actual", "harvest_simulated"
        elseif stuff == "emergence"
            "emergence_actual", "emergence_simulated"
        else
            error("unknown stuff: $stuff")
        end
    end

    # Stable styles per station (3 known; cycles if more)
    uniq_stations = collect(stations)
    colors = [:steelblue, :orange, :forestgreen]
    markers = [:circle, :diamond, :utriangle]
    style = Dict(s => (colors[mod1(i, length(colors))], markers[mod1(i, length(markers))])
                 for (i, s) in enumerate(uniq_stations))

    f = Makie.Figure(size=(1200, 800))
    axes = Axis[]   # keep to build one pooled legend later
    ncols = 2

    for (i, stuff) in enumerate(stuffs)
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1

        ax = Makie.Axis(f[row, col];
            title="$(crop_type) simulated vs actual $stuff",
            xlabel="actual value",
            ylabel="simulated value"
        )
        push!(axes, ax)

        xmin = Inf
        xmax = -Inf

        x_s, y_s = select_stuff(stuff)

        for s in uniq_stations
            # skip if keys not present
            if !haskey(ddd[s], crop_type)
                continue
            end
            d = ddd[s][crop_type]
            if !(haskey(d, x_s) && haskey(d, y_s))
                continue
            end

            xx = d[x_s]
            yy = d[y_s]
            x, y = CropGrowthTutorial.keep_only_notmissing(xx, yy)
            isempty(x) && continue

            xmin = min(xmin, minimum(x))
            xmax = max(xmax, maximum(x))

            c, m = style[s]
            Makie.scatter!(ax, x, y; color=c, marker=m, markersize=7, label=string(s))
        end

        # x = y reference line (no label → won’t show in pooled legend)
        if isfinite(xmin) && isfinite(xmax)
            Makie.lines!(ax, [xmin, xmax], [xmin, xmax]; color=:tomato, linestyle=:dash, label=nothing)
            Makie.xlims!(ax, xmin, xmax)
            Makie.ylims!(ax, xmin, xmax)
        end
    end

    # Single legend at the bottom spanning both columns, with one entry per station
    legend_elems = [MarkerElement(color=style[s][1], marker=style[s][2], markersize=10) for s in uniq_stations]
    legend_labels = string.(uniq_stations)
    f[3, 1:2] = Makie.Legend(
        f, legend_elems, legend_labels;
        orientation=:horizontal,   # → one row
        framevisible=false,
        halign=:center
    )

    return f
end
