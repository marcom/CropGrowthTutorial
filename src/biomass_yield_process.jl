"""
    calibrate_biomass_yield_parameters!(crop_dict, crop_name, soil_type, hk_clim_df, var_name,  phenology_df, target_array; kw...)

changes a canopy-root-precalibrated crop_dict with the calibrated crop parameters related to biomass (or yield) process (depending of var_name==:Biomass or :Yield) 
"""
function calibrate_biomass_yield_parameters!(crop_dict, crop_name, soil_type, hk_clim_df, var_name,  phenology_df, target_array; kw...)
    @assert var_name in [:Biomass, :Yield] "var_name should be either :Biomass or :Yield"
    
    # setup the crop since we need some parameters of it
    crop = AquaCrop.RepCrop()
    AquaCrop.set_crop!(crop, crop_name; aux = nothing)

    i = find_closest_to_median(target_array)

    crop_date = phenology_df[i,:sowingdate]
    target_value = target_array[i]

    u, pars = _calibrate_biomass_yield_parameters(crop, crop_dict, crop_name, soil_type, crop_date, hk_clim_df, target_value, var_name; kw...)
    
    for key in keys(pars)
        ii = Int(pars[key][1])
        crop_dict[key] = u[ii]
    end

    # given GDDCDC recalculate GDDaysToSenescence
    cdc = crop_dict["GDDCDC"]
    ltoharvest = crop_dict["GDDaysToHarvest"]
    ccx = crop.CCx
    ltosenscence = recalculate_GDDaysToSenescence(cdc, ltoharvest, ccx)
    crop_dict["GDDaysToSenescence"] = round(Int, ltosenscence)

    # reset HI to be an integer
    if var_name == :Yield
        crop_dict["HI"] = round(Int, crop_dict["HI"])
    end

    return nothing 
end

function _calibrate_biomass_yield_parameters(crop, crop_dict, crop_name, soil_type, crop_date, hk_clim_df, target_value, var_name; kw...)

    # parameters to calibrate have some bounds and initial condition picked by hand 
    pars = Dict(
        "GDtranspLow" => [1, 2, crop.Tupper-crop.Tbase-5, 2],
        "WP" => [2, 0.8*crop.WP, 1.33*crop.WP, 1.33*crop.WP],
        "GDDCGC" => [3, 0.8*crop_dict["GDDCGC"], 1.33*crop_dict["GDDCGC"], crop_dict["GDDCGC"]],
        "GDDCDC" => [4, 0.8*crop_dict["GDDCDC"], 1.33*crop_dict["GDDCDC"], crop_dict["GDDCDC"]*1.33],
    )
    if var_name == :Yield
        pars["HI"] = [5, 0.9*crop.HI, 1.1*crop.HI, crop.HI]
    end

    # optimization function
    function get_val(x, p)
        # given GDDCDC recalculate GDDaysToSenescence
        cdc = x[Int(p["GDDCDC"][1])]
        ltoharvest = crop_dict["GDDaysToHarvest"]
        ccx = crop.CCx
        ltosenscence = recalculate_GDDaysToSenescence(cdc, ltoharvest, ccx)
        _crop_dict = Dict(
                    "GDDaysToSenescence" => round(Int,ltosenscence),
                    "GDtranspLow" => x[Int(p["GDtranspLow"][1])],
                    "WP" => x[Int(p["WP"][1])],
                    "GDDCGC" => x[Int(p["GDDCGC"][1])],
                    "GDDCDC" => x[Int(p["GDDCDC"][1])],
        )
        if var_name == :Yield
            _crop_dict["HI"] = round(Int,x[Int(p["HI"][1])])
        end

        kw = (
            crop_dict = merge(crop_dict, _crop_dict),
            )


        cropfield, all_ok = run_simulation(soil_type, crop_name, crop_date, hk_clim_df; kw...)
        if !all_ok.logi
            return Inf
        else
            if var_name == :Yield
                return (cropfield.dayout[end,"Y(fresh)"].val - target_value)^2
            elseif var_name == :Biomass
                return (cropfield.dayout[end,"Biomass"].val - target_value)^2
            end
        end
    end

    # make the bounds
    lb = zeros(length(pars))
    ub = zeros(length(pars))
    u0 = zeros(length(pars))
    for key in keys(pars)
        ii = Int(pars[key][1])
        lb[ii] = pars[key][2]
        ub[ii] = pars[key][3]
        u0[ii] = (ub[ii] + lb[ii])/2 
    end

    # make optimization problem
    prob = OptimizationProblem(get_val, u0, pars, lb = lb, ub = ub)

    # solve the optimization problem
    maxiters = get(kw, :maxiters, 100)
    rng_seed = get(kw, :rng_seed, nothing)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited();
                maxiters=maxiters,
                (isnothing(rng_seed) ? (; ) : (; RngSeed=rng_seed, RandomizeRngSeed=false))...
    )
    return sol.u, pars
end

function recalculate_GDDaysToSenescence(cdcval, ltoharvest, ccxval)
    return ltoharvest - ((ccxval+2.29)/((cdcval)*3.33))*log(1 + 1/0.05) 
end

function find_closest_to_median(target_array)
    mu = mean(target_array)
    m = median(target_array)

    if length(target_array)%2 != 0
        _, i = findmin(xx->abs(xx-m), target_array) 
    else
        if mu >= m
            ii = -1
        else
            ii = 1
        end
        _, i = findmin(xx->abs(xx-m*(1 + ii*0.001) ), target_array) 
    end

    return i 
end
