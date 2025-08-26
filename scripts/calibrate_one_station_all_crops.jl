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
# using Infiltrator

Plots.gr()
Plots.default(
         fontfamily = "helvetica",
         guidefont  = font(12),   # axes titles
         tickfont   = font(10),   # tick labels
         legendfont = font(10),    # legend text
         guidefontfamily = "helvetica",
         fontfamily_subplot = "helvetica",
         legend_font_family = "helvetica",
         legend_title_font_family = "helvetica",
         plot_titlefontfamily = "helvetica",
         titlefontfamily = "helvetica",
         xguidefontfamily = "helvetica",
         yguidefontfamily = "helvetica",
         xtickfontfamilya = "helvetica",
         ytickfontfamilya = "helvetica",
         )


########
# MAIN #
########
function main(station_=nothing)
    # DFAULT STATION #
    default_station = "Jena"
    plottype = CropGrowthTutorial.PlotsPlotOption()

    ## READ DATA ##

    # Stations information
    stations_df = CropGrowthTutorial.get_phenology_stations()
    # Crop type general information
    general_crop_data_df = CropGrowthTutorial.get_crop_parameters()


    if isnothing(station_)
        println("\n\nWARNING: no station name given, selecting default station ",default_station)
        station_ = default_station
    elseif !(station_ in stations_df.station_name) 
        println("\n\nWARNING: station name ",station_," not found in satations_df, selecting default station ",default_station)
        station_ = default_station 
    end


    # find given station
    row = filter(x->x.station_name==station_, stations_df)[1,:]
    station_name = row.station_name
    soil_type = ismissing(row.soil_type) ? "missing" : row.soil_type
    fit_type = :unknown
    println("\n\nCCCCCCCCCCCCCCCCCCCCCCC\n","station_name ",station_name)

    # get climate data for given station
    hk_clim_df =  CropGrowthTutorial.get_climate_data(station_name)
    if isnothing(hk_clim_df)
        println("WARNING: No climate data for station ",station_name," going to next station")
        return nothing
    end
    # get yield data for given station
    hk_yield_df = CropGrowthTutorial.get_yield_data(station_name)
    if isnothing(hk_yield_df)
        println("WARNING: No yield data for station ",station_name," we will only make phenolgy and canopy fit")
        fit_type = :PhenologyCanopy 
    end

    # for each crop
    for roww in eachrow(general_crop_data_df)
        crop_type = roww.crop_type
        crop_name = roww.aquacrop_name
        plantingdens = roww.plantingdens
        println("\n###\n","crop_type ",crop_type)

        # Crop's phenology raw data for a given station
        phenology_raw_df = CropGrowthTutorial.get_crop_phenology_data(crop_type, station_name)
        if isnothing(phenology_raw_df)
            println("WARNING: No phenology data for station ",station_name,
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
        vv = unique(phenology_raw_df[!,:phase_id])
        if (!(sowing_phase in vv) || !(harvest_phase in vv))
            println("WARNING: No enough phenology data for station ",station_name,
                        " in crop_type ", crop_type,
                        " going to next crop")
            continue
        end


        # create a kw tuple with the additional information that we wish to pass to AquaCrop
        kw = ( crop_dict = Dict{String,Any}("PlantingDens" => plantingdens),);


        ## CALIBRATION ##
        if fit_type == :unknown
            println("WARNING: Unknown fit type selected for station ",station_name,
                        " in crop_type ", crop_type,
                        " going to next crop")
            continue
        else
            # 1. Calibrate Phenology 
            crop_dict, phenology_df = CropGrowthTutorial.calibrate_phenology_parameters(phenology_raw_df, crop_name, hk_clim_df, 
                                                              sowing_phase, harvest_phase; kw...)
            if isnothing(phenology_df)
                println("WARNING: No enough phenology data for station ",station_name,
                            " in crop_type ", crop_type,
                            " going to next crop")
                continue
            end
            

            # 2. Calibrate the canopy growth parameters
            CropGrowthTutorial.calibrate_canopy_root_parameters!(crop_dict, crop_name);

            if fit_type == :PhenologyCanopy
                fit_description =  "phenology canopy calibration"
            else
                # gather the yield/biomass data
                # NOTE: we assume we have this data, otherwise we should check and if missing go back to :PhenologyCanopy fit
                yield_df = filter(row->row.crop_type==crop_type, hk_yield_df);
                target_output = []
                if size(yield_df,1) == 0
                    println("WARNING: No yield data for crop_type="*crop_type*"on station_name="*
                        station_name*" going to next crop")
                    continue
                end

                for rowi in eachrow(phenology_df)
                    yeari = rowi.year
                    append!(target_output, yield_df[1,string(yeari)]/10)
                end

                # 3. Calibrate biomass-yield parameters
                CropGrowthTutorial.calibrate_biomass_yield_parameters!(crop_dict, crop_name, soil_type, hk_clim_df, fit_type,  phenology_df, target_output)

                if fit_type == :Biomass
                    fit_description =  "phenology canopy biomass calibration"
                elseif fit_type == :Yield
                    fit_description =  "phenology canopy yield calibration"
                end
            end

            ## PLOT RESULTS ##
            # plot of days to harvest actual and simulated
            f = CropGrowthTutorial.plot_correlation(plottype, phenology_df.harvest_actualdays, phenology_df.harvest_simulateddays, crop_type, station_name, "days to harvest")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages",crop_type*"_"*station_name*"harvest_correlation.svg"))
            else
                Makie.save(plotsdir(crop_type*"_"*station_name*"_correlation.png"), f)
            end
            f = CropGrowthTutorial.plot_correlation(plottype, phenology_df.emergence_actualdays, phenology_df.emergence_simulateddays, crop_type, station_name, "days to emergence")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages",crop_type*"_"*station_name*"emergence_correlation.svg"))
            else
                Makie.save(plotsdir(crop_type*"_"*station_name*"_correlation.png"), f)
            end
            f = CropGrowthTutorial.plot_correlation(plottype, phenology_df.beginflowering_actualdays, phenology_df.beginflowering_simulateddays, crop_type, station_name, "days to begin flowering")
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages",crop_type*"_"*station_name*"beginflowering_correlation.svg"))
            else
                Makie.save(plotsdir(crop_type*"_"*station_name*"_correlation.png"), f)
            end
            
            # f = CropGrowthTutorial.plot_GDD_stats_years(phenology_df, station_name)
            # Makie.save(plotsdir(crop_type*"_"*station_name*"_phenology.png"), f)

            # plot of crop growth simulation 
            i = CropGrowthTutorial.find_closest_to_median(target_output)
            sowingdate = phenology_df[i,"sowingdate"]
            kw = (
                crop_dict = crop_dict,
            )
            cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df; kw...)
            if all_ok.logi == false
                println("\nWARNING: bad simulation result on ", sowingdate)
                println("  ",all_ok.msg)
                continue
            end
            f = CropGrowthTutorial.plot_daily_stuff_one_year(plottype, cropfield, crop_type, soil_type; kw...)
            if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                Plots.savefig(f, plotsdir("svgimages", crop_type*"_"*station_name*"_cropevolution.svg"))
            else
                Makie.save(plotsdir(crop_type*"_"*station_name*"_cropevolution.png"), f)
            end

            # plot yield results
            if fit_type in [:Biomass, :Yield]
                years = phenology_df[:,"year"]
                simulated_yield = []
                simulated_yield_baseline = []
                for sowingdate in phenology_df[:,"sowingdate"]
                    # fited yield
                    cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df; kw...)
                    if all_ok.logi == false
                        println("\nWARNING: bad simulation result on ",sowingdate," for fited yield")
                        println("  ",all_ok.msg)
                        append!(simulated_yield, missing)
                        continue
                    end
                    if fit_type == :Biomass
                        append!(simulated_yield, ustrip(cropfield.dayout[end,"Biomass"]))
                    elseif fit_type == :Yield
                        append!(simulated_yield, ustrip(cropfield.dayout[end,"Y(fresh)"]))
                    end

                    # base line yield
                    cropfield, all_ok = CropGrowthTutorial.run_simulation(soil_type, crop_name, sowingdate, hk_clim_df)
                    if all_ok.logi == false
                        println("\nWARNING: bad simulation result on ",sowingdate," for baseline yield")
                        println("  ",all_ok.msg)
                        append!(simulated_yield_baseline, missing)
                        continue
                    end
                    if fit_type == :Biomass
                        append!(simulated_yield_baseline, ustrip(cropfield.dayout[end,"Biomass"]))
                    elseif fit_type == :Yield
                        append!(simulated_yield_baseline, ustrip(cropfield.dayout[end,"Y(fresh)"]))
                    end
                end

                dats = [ target_output, simulated_yield, simulated_yield_baseline ]
                cols = [ "actual yield", "fitted yield", "baseline yield" ]

                f = CropGrowthTutorial.plot_correlation(plottype, target_output, simulated_yield, crop_type, station_name, "yield data")
                if typeof(plottype) == CropGrowthTutorial.PlotsPlotOption
                    Plots.savefig(f, plotsdir("svgimages",crop_type*"_"*station_name*"yield_correlation.svg"))
                else
                    Makie.save(plotsdir(crop_type*"_"*station_name*"_correlation.png"), f)
                end

                f = CropGrowthTutorial.plot_yearly_data(dats, cols, years, crop_type, soil_type, station_name)
                Makie.save(plotsdir(crop_type*"_"*station_name*"_yearlydata.png"), f)
            end
            println("INFO: saved plots of crop_type "*crop_type*" from station_name "*station_name)
        end
    end
end

function write_fited_crop(crop_type, soil_type, crop_name, crop_dict, station_name, fit_description)
    # file where to save result
    file_name = crop_type*"_"*station_name*".TOML"
    full_filename = datadir("sims",file_name)

    # make crop
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux = crop_dict)

    # make header
    head = "Crop fited: crop_type="*crop_type*" station_name="*station_name*
            " soil_type="*soil_type*" "*fit_description

    # save crop
    AquaCrop.save_crop(full_filename, crop, head)
    return nothing
end
