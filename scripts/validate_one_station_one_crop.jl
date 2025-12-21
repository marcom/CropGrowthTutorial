###############
# DESCRIPTION #
###############
"""
    Validation script: simulate with calibrated TOML files

    Given a path file containing calibrated crop TOML and a station name,
    this script:
    1. Loads station/climate/yield data
    1. Extracts sowing dates per year from phenology_raw_df
    1. Runs simulations using the calibrated TOML file
    1. Compares simulated results vs actual yield (plot_correlation)
    1. Produces yearly plots (plot_yearly_data)
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
function main(crop_file_path, station_, _crop_type)
    plottype = CropGrowthTutorial.MakiePlotOption()

    ## READ DATA ##
    stations_df = CropGrowthTutorial.get_phenology_stations()
    general_crop_data_df = CropGrowthTutorial.get_crop_parameters()

    if !(station_ in stations_df.station_name)
        println("\n\nWARNING: station name ", station_, " not found in stations_df")
        return nothing
    end

    # find given station
    row = filter(x -> x.station_name == station_, stations_df)[1, :]
    station_name = row.station_name
    soil_type = ismissing(row.soil_type) ? "missing" : row.soil_type
    println("\n\n[VALIDATION] station_name ", station_name)

    # climate and yield
    hk_clim_df = CropGrowthTutorial.get_climate_data(station_name)
    if isnothing(hk_clim_df)
        println("WARNING: No climate data for station ", station_name, " aborting validation")
        return nothing
    end
    hk_yield_df = CropGrowthTutorial.get_yield_data(station_name)
    if isnothing(hk_yield_df)
        println("WARNING: No yield data for station ", station_name, " aborting validation")
        return nothing
    end

    # for each crop
    for roww in eachrow(general_crop_data_df)
        crop_type = roww.crop_type
        # skips rows that are not the wanted crop
        if crop_type != _crop_type
            continue
        else
            crop_name = roww.aquacrop_name
            plantingdens = roww.plantingdens
            println("\n### VALIDATING crop_type ", crop_type)

            # phenology raw data for the station/crop
            phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
            if isnothing(phenology_raw_df)
                println("WARNING: No phenology data for station ", station_name,
                    " in crop_type ", crop_type,
                    " skipping")
                continue
            end

            # Determine phenology phases (as in calibration)
            sowing_phase = 10
            harvest_phase = startswith(crop_type, "silage") ? 39 : 24
            vv = unique(phenology_raw_df[!, :phase_id])
            if (!(sowing_phase in vv) || !(harvest_phase in vv))
                println("WARNING: Not enough phenology data for station ", station_name,
                    " in crop_type ", crop_type, " skipping")
                continue
            end

            # Yield row for this crop (wide format: one column per year)
            yield_df = filter(r -> r.crop_type == crop_type, hk_yield_df)
            if size(yield_df, 1) == 0
                println("WARNING: No yield data for crop_type=", crop_type, " on station_name=", station_name, " skipping")
                continue
            end

            # Path to calibrated TOML for this crop/station
            if !isfile(crop_file_path)
                println("WARNING: Calibrated TOML file not found: ", crop_file_path, " skipping crop")
                continue
            end

            # create empty crop and read the file
            c = AquaCrop.RepCrop()
            AquaCrop.load_gvars_from_toml!(c, crop_file_path)

            # Convert RepCrop struct to a Dict and force planting density from general data
            acceptd_fields = [:GDDCGC, :GDDaysToGermination, :GDDaysToFlowering, :GDDaysToMaxRooting,
                :GDDLengthFlowering, :GDDCDC, :RootMax, :GDDaysToHarvest, :GDDaysToFullCanopy,
                :GDtranspLow, :PlantingDens, :GDDaysToSenescence, :WP]
            crop_dict = Dict{String,Any}(string(k) => getfield(c, k) for k in fieldnames(typeof(c)) if k in acceptd_fields)
            crop_dict["PlantingDens"] = plantingdens

            # kw: pass the full crop dict
            kw_sim = (crop_dict=crop_dict,)

            # Build phenology_df using calibration phenology step (we only use the dataframe)
            # Use PlantingDens as in calibration kw
            pheno_kw = (crop_dict=Dict{String,Any}("PlantingDens" => plantingdens),)
            _crop_dict_ignore, phenology_df = CropGrowthTutorial.calibrate_phenology_parameters(
                phenology_raw_df, crop_name, hk_clim_df, sowing_phase, harvest_phase; pheno_kw...
            )
            if isnothing(phenology_df)
                println("WARNING: Not enough processed phenology for station ", station_name,
                    " in crop_type ", crop_type, " skipping")
                continue
            end

            # Collect target outputs and matched sowing dates (use processed phenology_df; see calibration’s 169–198)
            target_output = Float64[]
            _phenology_df = copy(phenology_df)
            inds = Int[]

            for (_i, rowi) in enumerate(eachrow(phenology_df))
                yeari = rowi.year
                crop_date = rowi.sowingdate
                if (hk_clim_df[1, "date"] < crop_date < hk_clim_df[end, "date"]) &&
                   (string(yeari) in names(yield_df)) &&
                   !(ismissing(yield_df[1, string(yeari)]))
                    append!(target_output, yield_df[1, string(yeari)] / 10)
                else
                    push!(inds, _i)
                end
            end
            if length(inds) > 0
                deleteat!(_phenology_df, inds)
            end

            if length(target_output) == 0
                println("WARNING: Not enough yield data after matching years/sowing for crop_type=", crop_type, " on station_name=", station_name, " skipping")
                continue
            end

            # Phenology correlation plots (harvest, emergence, beginflowering) as in calibration
            f_h = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.harvest_actualdays, _phenology_df.harvest_simulateddays, crop_type, station_name, "days to harvest")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f_h, plotsdir("svgimages", crop_type * "_" * station_name * "harvest_correlation_validation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "harvest_correlation_validation.svg"), f_h)
                Makie.save(plotsdir(crop_type * "_" * station_name * "harvest_correlation_validation.png"), f_h)
            end

            f_e = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.emergence_actualdays, _phenology_df.emergence_simulateddays, crop_type, station_name, "days to emergence")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f_e, plotsdir("svgimages", crop_type * "_" * station_name * "emergence_correlation_validation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "emergence_correlation_validation.svg"), f_e)
                Makie.save(plotsdir(crop_type * "_" * station_name * "emergence_correlation_validation.png"), f_e)
            end

            f_b = CropGrowthTutorial.plot_correlation(plottype, _phenology_df.beginflowering_actualdays, _phenology_df.beginflowering_simulateddays, crop_type, station_name, "days to begin flowering")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f_b, plotsdir("svgimages", crop_type * "_" * station_name * "beginflowering_correlation_validation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "beginflowering_correlation_validation.svg"), f_b)
                Makie.save(plotsdir(crop_type * "_" * station_name * "beginflowering_correlation_validation.png"), f_b)
            end

            # Simulate for each sowing date with calibrated TOML
            simulated_yield = Vector{Union{Missing,Float64}}()
            use_biomass = startswith(crop_type, "silage")
            for sowingdate in _phenology_df[:, "sowingdate"]
                cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df; kw_sim...)
                if all_ok.logi == false
                    println("\nWARNING: bad simulation result on ", sowingdate)
                    println("  ", all_ok.msg)
                    push!(simulated_yield, missing)
                    continue
                end
                if use_biomass
                    push!(simulated_yield, ustrip(cropfield.dayout[end, "Biomass"]))
                else
                    push!(simulated_yield, ustrip(cropfield.dayout[end, "Y(fresh)"]))
                end
            end

            # Plots: correlation and yearly series
            f1 = CropGrowthTutorial.plot_correlation(plottype, target_output, simulated_yield, crop_type, station_name, "yield data")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f1, plotsdir("svgimages", crop_type * "_" * station_name * "yield_correlation_validation.svg"))
            else
                Makie.save(plotsdir("svgimages", crop_type * "_" * station_name * "yield_correlation_validation.svg"), f1)
                Makie.save(plotsdir(crop_type * "_" * station_name * "yield_correlation_validation.png"), f1)
            end

            dats = [target_output, simulated_yield]
            cols = ["actual yield", "simulated yield"]
            years = _phenology_df[:, "year"]
            f2 = CropGrowthTutorial.plot_yearly_data(dats, cols, years, crop_type, soil_type, station_name)
            Makie.save(plotsdir(crop_type * "_" * station_name * "_yearlydata_validation.png"), f2)

            println("INFO: saved validation plots of crop_type ", crop_type, " from station_name ", station_name)
        end

    end
    return nothing
end
