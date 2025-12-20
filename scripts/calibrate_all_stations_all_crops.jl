###############
# DESCRIPTION #
###############
"""
    We make the crop calibration for all the crops and all the stations

    Note:
    1. Some stations do not have enough phenology data to calibrate the phenology
    1. Some stations do not have yield/biomass data

    If we only have full data for a given crop in a given station we do the full calibration
    If we only have enough data for the phenology calibration we only do this
    If we do not have enough data for anything we skip the crop

    Results of crop calibration are on `data/sims`
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
# using Infiltrator


########
# MAIN #
########
function main()

    ## READ DATA ##

    # Stations information
    stations_df = CropGrowthTutorial.get_phenology_stations()
    # Crop type general information
    general_crop_data_df = CropGrowthTutorial.get_crop_parameters()

    # for each station
    for row in eachrow(stations_df)
        station_name = row.station_name
        soil_type = ismissing(row.soil_type) ? "missing" : row.soil_type
        fit_type = :unknown
        println("\n\nCCCCCCCCCCCCCCCCCCCCCCC\n", "station_name ", station_name)

        # get climate data for given station
        hk_clim_df = CropGrowthTutorial.get_climate_data(station_name)
        if isnothing(hk_clim_df)
            println("WARNING: No climate data for station ", station_name, " going to next station")
            continue
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

                ## WRITE RESULTS ##
                write_fited_crop(crop_type, soil_type, crop_name, crop_dict, station_name, fit_description)
                println("INFO: saved crop_type " * crop_type * " from station_name " * station_name)
            end
        end
    end
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
