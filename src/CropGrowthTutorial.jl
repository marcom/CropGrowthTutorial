module CropGrowthTutorial 
# ---- imports ----
using DrWatson
using Printf
using StatsBase
using DataFrames
using Dates
using CSV
using AquaCrop
using CairoMakie
using Plots
using Unitful
using Optimization
using OptimizationBBO
using Infiltrator

# ---- includes ----
include("global_data.jl")
include("read_data.jl")
include("plot_data.jl")
include("plot_data_plots.jl")
include("simulation.jl")
include("phenology_process.jl")
include("canopy_root_process.jl")
include("biomass_yield_process.jl")


#---- exports ----
export   head
#        get_climate_data,
#        list_climate_files,
#        get_yield_data,
#        list_yield_data_files,
#        get_crop_phenology_data,
#        list_phenology_data_files,
#        get_phenology_stations,
#        get_all_phenology_stations,
#        get_all_phenology_phases,
#        get_crop_parameters,
#        process_crop_phenology,
#        run_simulation,
#        run_cropgdd,
#        process_crop_phenology_gdd

end

