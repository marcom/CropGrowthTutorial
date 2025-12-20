##########################
# functions to plot data #
##########################

# To handle multiple dispatch
abstract type AbstractPlotOption end

struct MakiePlotOption <: AbstractPlotOption end
# this is de default option

struct PlotsPlotOption <: AbstractPlotOption end
# note that we only have Plots recipies for plot_daily_stuff_one_year and plot_correlation


"""
    plot_daily_stuff_one_year(cropfield::AquaCropField, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass" ], plot_label=true; kw...)

plots the data in `cropfield.dayout` considering the columns `cols`
"""
function plot_daily_stuff_one_year(plottype::T, cropfield::AquaCropField, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass"], plot_label=true; kw...) where {T<:AbstractPlotOption}
    plot_daily_stuff_one_year(plottype, cropfield.dayout, crop_type, soil_type, cols, plot_label; kw...)
end

function plot_daily_stuff_one_year(cropfield::AquaCropField, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass"], plot_label=true; kw...)
    plot_daily_stuff_one_year(MakiePlotOption(), cropfield.dayout, crop_type, soil_type, cols, plot_label; kw...)
end

function plot_daily_stuff_one_year(cropfield::AbstractDataFrame, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass"], plot_label=true; kw...)
    plot_daily_stuff_one_year(MakiePlotOption(), cropfield, crop_type, soil_type, cols, plot_label; kw...)
end

function plot_daily_stuff_one_year(plottype::MakiePlotOption, cropfield::AbstractDataFrame, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass"], plot_label=true; kw...)
    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end


    x = cropfield[1:xx, "Date"]
    aux_sz = round(Int, sqrt(length(cols)))

    f = Makie.Figure()
    for (i, coli) in enumerate(cols)
        ii, jj = divrem(i - 1, aux_sz)
        ax = Makie.Axis(f[jj, ii],
            title=coli * " vs Date",
            xlabel="Date",
            ylabel=coli,
        )
        Makie.lines!(ax, x, ustrip.(cropfield[1:xx, coli]))

        # to plot a vertical line on special days
        if haskey(kw, :end_day)
            Makie.vlines!(ax, Makie.date_to_number(DateTime, kw[:end_day]), color=:tomato, linestyle=:dash, label="harvest day")
        end
        if haskey(kw, :emergence_day) && !ismissing(kw[:emergence_day])
            Makie.vlines!(ax, Makie.date_to_number(DateTime, kw[:emergence_day]), color=:green, linestyle=:dash, label="emergence day")
        end
        if haskey(kw, :beginflowering_day) && !ismissing(kw[:beginflowering_day])
            Makie.vlines!(ax, Makie.date_to_number(DateTime, kw[:beginflowering_day]), color=:yellow, linestyle=:dash, label="flowering day")
        end
        if haskey(kw, :endflowering_day) && !ismissing(kw[:endflowering_day])
            Makie.vlines!(ax, Makie.date_to_number(DateTime, kw[:endflowering_day]), color=:orange, linestyle=:dash, label="end flowering day")
        end

        # to plot a horizontal line in a minimum temperature
        if (coli == "Tmin") & haskey(kw, :Tmin)
            Makie.hlines!(ax, kw[:Tmin], color=:tomato, linestyle=:dash, label="min T flowering")
        end
        # Makie.axislegend(position = :rb)

        ax.xticklabelrotation = π / 4
        ax.xticklabelsize = 8
        ax.yticklabelsize = 8
    end

    if plot_label
        # Add a general title
        Makie.Label(f[-1, :], "Simulation results for " * crop_type, fontsize=15, valign=:center, font=:bold)

        # Add small text (e.g., as a footer)
        Makie.Label(f[aux_sz, :], "Simulation starts in " * string(cropfield[1, "Date"]) * " and ends in " * string(cropfield[end, "Date"]) * " soil_type is " * soil_type, fontsize=12, valign=:bottom)
    end

    return f
end

"""
    plot_yearly_data(datas, cols, years_interval, crop_type, soil_type, region_name)

Plots the results in `datas::Vector{Vector{Float64}}` vs the years in `years_interval`.
Uses the names in `cols::Vector{String}` for the labels
"""
function plot_yearly_data(datas, cols, years_interval, crop_type, soil_type, region_name)
    @assert length(datas) == length(cols) "Bad number of data or cols names"
    x = years_interval

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title=crop_type * " yield over time in region " * region_name,
        xlabel="Year",
        ylabel="ton/ha"
    )
    for (i, coli) in enumerate(cols)
        Makie.scatterlines!(ax, x, datas[i], label=coli)
    end
    f[1, 2] = Makie.Legend(f, ax, "soil_type is " * soil_type, titlefont=:regular, framevisible=false)

    return f
end

"""
    plot_gddaily_stuff_one_year(cropfield::AquaCropField, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass" ])

plots the data in `cropfield.dayout` considering the columns `cols`, it uses Growing Degree Days on the `x` axis.
"""
function plot_gddaily_stuff_one_year(cropfield::AquaCropField, crop_type, soil_type, cols=["CC", "Stage", "Y(fresh)", "Biomass"])
    x = cumsum(cropfield.dayout.GD)
    aux_sz = round(Int, sqrt(length(cols)))

    f = Makie.Figure()
    for (i, coli) in enumerate(cols)
        ii, jj = divrem(i - 1, aux_sz)
        ax = Makie.Axis(f[ii, jj],
            title=coli * " vs GDDay",
            xlabel="GDDay",
            ylabel=coli,
        )
        Makie.lines!(ax, x, ustrip.(cropfield.dayout[!, coli]))
    end
    # Add a general title
    Makie.Label(f[-1, :], "Simulation results for " * crop_type, fontsize=15, valign=:center, font=:bold)

    # Add small text (e.g., as a footer)
    Makie.Label(f[aux_sz, :], "Simulation starts in " * string(cropfield.dayout[1, "Date"]) * " and ends in " * string(cropfield.dayout[end, "Date"]) * " soil_type is " * soil_type, fontsize=12, valign=:bottom)

    return f
end

"""
    plot_correlation(x, y, crop_type, region_name, variable_name)

plots the simulated yield vs measured yield, it also finds the Pearson and Spearman correlation between them
"""
function plot_correlation(xx, yy, crop_type, region_name, variable_name)
    plot_correlation(MakiePlotOption(), xx, yy, crop_type, region_name, variable_name)
end

function plot_correlation(plottype::MakiePlotOption, xx, yy, crop_type, region_name, variable_name)
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title=crop_type * " simulated vs actual " * variable_name * " in region " * region_name,
        xlabel="actual value",
        ylabel="simulated value"
    )

    x, y = keep_only_notmissing(xx, yy)
    x = Float64.(x)
    y = Float64.(y)
    per = cor(x, y)
    spe = corspearman(x, y)
    mse = msd(x, y)
    x_bar = mean(x)
    willmott = 1 - mse * length(x) / mapreduce(x -> x^2, +, (abs.(x .- x_bar) + abs.(y .- x_bar)))
    Makie.scatter!(ax, x, y; label="data")
    sort!(x)
    Makie.lines!(ax, x, x; color=:tomato, linestyle=:dash, label="x=y")
    # ss = @sprintf "Pearson is %4.2f \nSpearman is %4.2f \n Willmott is %5.2f" per spe willmott
    # f[1, 2] = Makie.Legend(f, ax, ss, titlefont = :regular, framevisible = false)

    return f
end

"""
    plot_soil_wc(cropfield::AquaCropField; kw...)

plots the daily soil water content
"""
function plot_soil_wc(cropfield::AquaCropField; kw...)
    plot_soil_wc(cropfield.dayout; kw...)
end

function plot_soil_wc(cropfield::AbstractDataFrame; kw...)
    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end

    x = cropfield[1:xx, "Date"]

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title="Soil Water Content",
        xlabel="Date",
        ylabel="mm"
    )

    Makie.lines!(ax, x, ustrip.(cropfield[1:xx, "Wr"] - cropfield[1:xx, "Wr(SAT)"]), label="Water in root zone")
    Makie.lines!(ax, x, ustrip.(cropfield[1:xx, "Wr(FC)"] - cropfield[1:xx, "Wr(SAT)"]), linestyle=:dash, label="Field Capacity")
    Makie.lines!(ax, x, ustrip.(cropfield[1:xx, "Wr(PWP)"] - cropfield[1:xx, "Wr(SAT)"]), linestyle=:dash, label="Wilting Point")
    # Makie.lines!(ax, x, ustrip.(cropfield[1:xx,"Wr(exp)"] - cropfield[1:xx,"Wr(SAT)"]), linestyle=:dash, label="exp")
    # Makie.lines!(ax, x, ustrip.(cropfield[1:xx,"Wr(sto)"] - cropfield[1:xx,"Wr(SAT)"]), linestyle=:dash, label="sto")
    # Makie.lines!(ax, x, ustrip.(cropfield[1:xx,"Wr(sen)"] - cropfield[1:xx,"Wr(SAT)"]), linestyle=:dash, label="sen")


    Makie.axislegend(position=:lt)

    ax.xticklabelrotation = π / 4
    ax.xticklabelsize = 8
    ax.yticklabelsize = 8

    ax2 = Makie.Axis(f[1, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabel="mm", ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)

    Makie.barplot!(ax2, x, ustrip.(cropfield[1:xx, "Rain"]), label="Rain", gap=0, color=(:gray85, 0.5), strokecolor=:black, strokewidth=1)
    ax2.xticklabelrotation = π / 4
    ax2.xticklabelsize = 8
    ax2.yticklabelsize = 8
    Makie.axislegend(position=:rt)


    ax3 = Makie.Axis(f[1, 1])
    Makie.hidespines!(ax3)
    Makie.hidexdecorations!(ax3)
    Makie.hideydecorations!(ax3)
    Makie.lines!(ax3, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)

    return f
end

"""
    plot_eto_tr(cropfield::AquaCropField; kw...)

plots the daily potential and real crop's transpiration and the FAO Evapotranspiration
"""
function plot_eto_tr(cropfield::AquaCropField; kw...)
    plot_eto_tr(cropfield.dayout; kw...)
end

function plot_eto_tr(cropfield::AbstractDataFrame; kw...)

    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end

    x = cropfield[1:xx, "Date"]

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title="Evapotranspiration",
        xlabel="Date",
        ylabel="mm"
    )
    Makie.barplot!(ax, x, ustrip.(cropfield[1:xx, "ETo"]), label="ETo", gap=0, color=:gray85, strokecolor=:black, strokewidth=1)
    Makie.barplot!(ax, x, ustrip.(cropfield[1:xx, "Trx"]), label="Potential Transpiration", gap=0, color=(:dodgerblue, 0.5), strokecolor=:black, strokewidth=1)
    Makie.barplot!(ax, x, ustrip.(cropfield[1:xx, "Tr"]), label="Actual Transpiration", gap=0, color=(:lightsalmon, 0.3), strokecolor=:black, strokewidth=1)
    ax.xticklabelrotation = π / 4
    ax.xticklabelsize = 8
    ax.yticklabelsize = 8
    Makie.axislegend(position=:lt)


    ax2 = Makie.Axis(f[1, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.hideydecorations!(ax2)
    # Makie.lines!(ax2, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)
    Makie.lines!(ax2, x, cropfield[1:xx, "CC"], color=:gold, linestyle=:dash)

    return f
end

"""
    plot_crop_stress(cropfield::AquaCropField; kw...)

plots the crop's daily stresses
"""
function plot_crop_stress(cropfield::AquaCropField, var::Symbol=:all; kw...)
    plot_crop_stress(cropfield.dayout, var; kw...)
end

function plot_crop_stress(cropfield::AbstractDataFrame, var::Symbol=:all; kw...)
    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end

    x = cropfield[1:xx, "Date"]

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title="Crop Stresses",
        xlabel="Date",
        ylabel="%"
    )
    if var == :all
        plt_wtr = true
        plt_tem = true
    elseif var == :water
        plt_wtr = true
        plt_tem = false
    elseif var == :temperature
        plt_wtr = false
        plt_tem = true
    end

    if plt_wtr
        Makie.lines!(ax, x, cropfield[1:xx, "StExp"], label="Expansion stress")
        Makie.lines!(ax, x, cropfield[1:xx, "StSto"], label="Stomatal stress")
        Makie.lines!(ax, x, cropfield[1:xx, "StSen"], label="Senescence stress")
        # Makie.lines!(ax, x, cropfield[1:xx, "StSalt"], label="Salt stress")
        # Makie.lines!(ax, x, cropfield[1:xx, "StWeed"], label="StWeed")
    end
    if plt_tem
        Makie.lines!(ax, x, cropfield[1:xx, "StTr"], label="Temperature stress")
    end
    Makie.axislegend(position=:lt)

    ax.xticklabelrotation = π / 4
    ax.xticklabelsize = 8
    ax.yticklabelsize = 8
    Makie.ylims!(ax, -1, nothing)

    ax2 = Makie.Axis(f[1, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.hideydecorations!(ax2)
    Makie.lines!(ax2, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)
    Makie.ylims!(ax2, 0, nothing)

    return f
end

"""
    plot_wpi(cropfield::AquaCropField; kw...)

plots the crop's daily water productivity
"""
function plot_wpi(cropfield::AquaCropField; kw...)
    plot_wpi(cropfield.dayout, cropfield.crop.WP, cropfield.crop.WPy; kw...)
end

function plot_wpi(cropfield::AbstractDataFrame, wp=100, wpy=100; kw...)
    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end

    x = cropfield[1:xx, "Date"]

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title="Water productivity",
        xlabel="Date",
        ylabel="%"
    )
    Makie.lines!(ax, x, ustrip.(cropfield[1:xx, "WP"]), label="Water Productivity")
    Makie.hlines!(ax, wp, color=:tomato, linestyle=:dash, label="target WP")
    Makie.axislegend(position=:rb)

    ax.xticklabelrotation = π / 4
    ax.xticklabelsize = 8
    ax.yticklabelsize = 8
    Makie.ylims!(ax, -1, nothing)

    ax2 = Makie.Axis(f[1, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.hideydecorations!(ax2)
    Makie.lines!(ax2, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)

    return f
end

"""
    plot_temp(cropfield::AquaCropField; kw...)

plots the daily tempertaure and the crop's temperature parameters
"""
function plot_temp(cropfield::AquaCropField; kw...)
    plot_temp(cropfield.dayout, cropfield.crop.Tcold, cropfield.crop.Theat, cropfield.crop.GDtranspLow; kw...)
end

function plot_temp(cropfield::AbstractDataFrame, tcold=10, theat=40, gd=12; kw...)
    # crop.Tcold = 10         # Minimum air temperature below which pollination starts to fail (cold stress) (°C)
    # crop.Theat = 40         # Maximum air temperature above which pollination starts to fail (heat stress) (°C)

    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end

    x = cropfield[1:xx, "Date"]

    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1],
        title="Temperature",
        xlabel="Date",
        ylabel="°C"
    )

    Makie.lines!(ax, x, ustrip.(uconvert.(u"°C", cropfield[1:xx, "Tavg"])), label="Tavg", color=:chartreuse3)
    Makie.lines!(ax, x, ustrip.(uconvert.(u"°C", cropfield[1:xx, "Tmin"])), color=:turquoise, linestyle=:dash)
    Makie.lines!(ax, x, ustrip.(uconvert.(u"°C", cropfield[1:xx, "Tmax"])), color=:tomato, linestyle=:dash)
    Makie.hlines!(ax, tcold, color=:turquoise2, linestyle=:dashdotdot, label="tcold flowering")
    Makie.hlines!(ax, theat, color=:tomato2, linestyle=:dashdotdot, label="theat flowering")

    Makie.axislegend(position=:lt)

    ax.xticklabelrotation = π / 4
    ax.xticklabelsize = 8
    ax.yticklabelsize = 8

    ax2 = Makie.Axis(f[1, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.hideydecorations!(ax2)
    Makie.lines!(ax2, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)


    ax = Makie.Axis(f[2, 1],
        title="GD",
        xlabel="Date",
        ylabel="°C"
    )
    Makie.lines!(ax, x, cropfield[1:xx, "GD"], label="GD", color=:chartreuse3)
    Makie.hlines!(ax, gd, color=:turquoise, linestyle=:dashdotdot, label="canopy growht limit")

    ax2 = Makie.Axis(f[2, 1], yticklabelcolor=:tomato, yaxisposition=:right, ylabelcolor=:tomato)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.hideydecorations!(ax2)
    Makie.lines!(ax2, x, cropfield[1:xx, "Stage"], color=:gold, linestyle=:dash)
    return f
end

"""
    plot_GDD_stats_violin(df::AbstractDataFrame, station_name)

plots the GDD statistics for different crop phenology phases and stations
"""
function plot_GDD_stats_violin(df::AbstractDataFrame, station_name)

    function plot_violin_col!(ax, col_1::AbstractArray)
        v_1 = collect(skipmissing(col_1))

        if length(v_1) > 0
            x_1 = fill(1, length(v_1))

            # Plot violins
            Makie.violin!(ax, x_1, v_1; show_median=true, datalimits=extrema)
        end
    end


    f = Makie.Figure()

    ax1 = Makie.Axis(f[1, 1],
        title="Harvest GDD",
        ylabel="GDD",
        yticklabelsize=8,
        xticks=([1], [station_name])
    )
    plot_violin_col!(ax1, df.harvest_actualgdd)
    ax1.xticklabelrotation = π / 4

    ax2 = Makie.Axis(f[2, 1],
        title="Begin Flowering GDD",
        ylabel="GDD",
        xticks=([1], [station_name])
    )
    plot_violin_col!(ax2, df.beginflowering_actualgdd)
    ax2.xticklabelrotation = π / 4

    ax3 = Makie.Axis(f[1, 2],
        title="End Flowering GDD",
        ylabel="GDD",
        xticks=([1], [station_name])
    )
    plot_violin_col!(ax3, df.endflowering_actualgdd)
    ax3.xticklabelrotation = π / 4

    ax4 = Makie.Axis(f[2, 2],
        title="Emergence GDD",
        ylabel="GDD",
        xticks=([1], [station_name])
    )
    plot_violin_col!(ax4, df.emergence_actualgdd)
    ax4.xticklabelrotation = π / 4
    return f
end

"""
    plot_GDD_stats_years(df::AbstractDataFrame, crop_type="")

plots the GDD and CD for different crop phenology phases and stations
"""
function plot_GDD_stats_years(df::AbstractDataFrame, station_name)
    f = Makie.Figure()

    ax1 = Makie.Axis(f[1, 1],
        title="Harvest Day difference",
        ylabel="Days",
        xlabel="Year"
    )
    Makie.scatter!(ax1, year.(df.sowingdate), (df.harvest_actualdays - df.harvest_simulateddays), label=station_name, color=1, alpha=0.8, colormap=:tab10, colorrange=(1, 10))

    ax1_1 = Makie.Axis(f[1, 1],
        ylabel="GDDays",
        yticklabelcolor=:tomato,
        yaxisposition=:right,
        ylabelcolor=:tomato,
    )
    Makie.scatter!(ax1_1, year.(df.sowingdate), (df.harvest_actualgdd - df.harvest_simulatedgdd), label=station_name, color=2, alpha=0.8, colormap=:tab10, colorrange=(1, 10), marker=:rtriangle)
    Makie.hidespines!(ax1_1)
    Makie.hidexdecorations!(ax1_1)

    ax2 = Makie.Axis(f[2, 1],
        title="Begin Flowering Day difference",
        ylabel="Days",
        xlabel="Year"
    )
    x, y = keep_only_notmissing(year.(df.sowingdate), (df.beginflowering_actualdays - df.beginflowering_simulateddays))
    if length(x) > 1
        Makie.scatter!(ax2, x, y, label=station_name, color=1, alpha=0.8, colormap=:tab10, colorrange=(1, 10))
    end
    ax2_1 = Makie.Axis(f[2, 1],
        ylabel="GDDays",
        yticklabelcolor=:tomato,
        yaxisposition=:right,
        ylabelcolor=:tomato,
    )
    x, y = keep_only_notmissing(year.(df.sowingdate), (df.beginflowering_actualgdd - df.beginflowering_simulatedgdd))
    if length(x) > 1
        Makie.scatter!(ax2_1, x, y, label="uhk", color=2, alpha=0.8, colormap=:tab10, colorrange=(1, 10), marker=:rtriangle)
    end
    Makie.hidespines!(ax2_1)
    Makie.hidexdecorations!(ax2_1)

    ax3 = Makie.Axis(f[3, 1],
        title="Emergence Day difference",
        ylabel="Days",
        xlabel="Year"
    )
    x, y = keep_only_notmissing(year.(df.sowingdate), (df.emergence_actualdays - df.emergence_simulateddays))
    if length(x) > 1
        Makie.scatter!(ax3, x, y, label="uhk", color=1, alpha=0.8, colormap=:tab10, colorrange=(1, 10))
    end
    ax3_1 = Makie.Axis(f[3, 1],
        ylabel="GDDays",
        yticklabelcolor=:tomato,
        yaxisposition=:right,
        ylabelcolor=:tomato,
    )
    x, y = keep_only_notmissing(year.(df.sowingdate), (df.emergence_actualgdd - df.emergence_simulatedgdd))
    if length(x) > 1
        Makie.scatter!(ax3_1, x, y, label="uhk", color=2, alpha=0.8, colormap=:tab10, colorrange=(1, 10), marker=:rtriangle)
    end
    Makie.hidespines!(ax3_1)
    Makie.hidexdecorations!(ax3_1)


    Makie.Legend(f[1, 2], ax1, "Days", titlefont=:regular, framevisible=false)
    Makie.Legend(f[2, 2], ax1_1, "GDDays", titlefont=:regular, framevisible=false, titlecolor=:tomato)

    return f
end

function keep_only_notmissing(x, y)
    @assert length(x) == length(y)
    xx = []
    yy = []
    for (x_, y_) in zip(x, y)
        if ismissing(x_) || ismissing(y_)
            continue
        else
            append!(xx, x_)
            append!(yy, y_)
        end
    end
    return xx, yy
end
