# Plots stuff
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


function plot_daily_stuff_one_year(
    plottype::PlotsPlotOption,
    cropfield::AbstractDataFrame,
    crop_type,
    soil_type,
    cols=["CC", "Stage", "Y(fresh)", "Biomass"],
    plot_label::Bool=true; kw...
)

    xx = findlast(x -> x == 4, cropfield.Stage) + 1
    d_ii = cropfield[xx, :Date]
    if haskey(kw, :end_day)
        d_ii_ = kw[:end_day] + Day(1)
        if d_ii_ > d_ii
            xx = findfirst(x -> x == d_ii_, cropfield.Date)
        end
    end


    x = cropfield[1:xx, "Date"]

    # Grid layout like your aux_sz ~= ceil(sqrt(n))
    n = length(cols)
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)

    p = Plots.plot(layout=(nrows, ncols), legend=false)


    for (i, coli) in enumerate(cols)
        y = ustrip(cropfield[1:xx, coli])

        Plots.plot!(p, x, y;
            subplot=i,
            xlabel="Date",
            ylabel=string(coli),
            title=string(coli, " vs Date"),
            lw=2,
            xrotation=45
        )

        # Vertical markers (Dates)
        if haskey(kw, :end_day)
            Plots.vline!(p, [kw[:end_day].value]; subplot=i, color=:tomato, linestyle=:dash, label="harvest day")
        end
        if haskey(kw, :emergence_day) && !ismissing(kw[:emergence_day])
            Plots.vline!(p, [kw[:emergence_day].value]; subplot=i, color=:green, linestyle=:dash, label="emergence day")
        end
        if haskey(kw, :beginflowering_day) && !ismissing(kw[:beginflowering_day])
            Plots.vline!(p, [kw[:beginflowering_day].value]; subplot=i, color=:yellow, linestyle=:dash, label="flowering day")
        end
        if haskey(kw, :endflowering_day) && !ismissing(kw[:endflowering_day])
            Plots.vline!(p, [kw[:endflowering_day].value]; subplot=i, color=:orange, linestyle=:dash, label="end flowering day")
        end

        # Horizontal marker for Tmin panel
        if string(coli) == "Tmin" && haskey(kw, :Tmin)
            Plots.hline!(p, [kw[:Tmin]]; subplot=i, color=:tomato, linestyle=:dash, label="min T flowering")
        end
    end

    if plot_label
        # Overall title for the full layout (keeps text editable in SVG)
        Plots.plot!(p; plot_title="Simulation results for $(crop_type)")
        # txt = "Simulation starts on $(date_col[1]) and ends on $(date_col[end]); soil_type: $(soil_type)"
        # Plots.annotate!(p, subplot = n; x[1], minimum(_ustrip(cropfield[1:xx, cols[end]])), text(txt, 8, :left))
    end

    return p
end

function plot_correlation(plottype::PlotsPlotOption, xx, yy, crop_type, region_name, variable_name)
    x, y = keep_only_notmissing(xx, yy)
    x = Float64.(x)
    y = Float64.(y)
    per = cor(x, y)
    spe = corspearman(x, y)
    mse = msd(x, y)
    x_bar = mean(x)
    willmott = 1 - mse * length(x) / mapreduce(x -> x^2, +, (abs.(x .- x_bar) + abs.(y .- x_bar)))

    # Layout: left = scatter/identity; right = text panel
    layout = Plots.@layout [A{0.74w} B]
    p = Plots.plot(layout=layout)

    # Left panel: scatter + x=y
    ttl = string(crop_type, " simulated vs actual ", variable_name, " in region ", region_name)
    minxy = min(minimum(x), minimum(y))
    maxxy = max(maximum(x), maximum(y))

    Plots.scatter!(p, x, y;
        subplot=1,
        label="data",
        xlabel="actual value",
        ylabel="simulated value",
        title=ttl
    )

    Plots.plot!(p, [minxy, maxxy], [minxy, maxxy];
        subplot=1,
        color=:tomato,
        linestyle=:dash,
        label="x = y"
    )

    Plots.xlims!(p, minxy, maxxy; subplot=1)
    Plots.ylims!(p, minxy, maxxy; subplot=1)

    # Right panel: multiline stats text (kept as <text> in SVG)
    # stats = @sprintf("Pearson:  %4.2f\nSpearman: %4.2f\nWillmott: %5.2f", per, spe, willmott)
    # Plots.plot!(p; subplot=2, framestyle=:none, xticks=false, yticks=false, xlim=(0,1), ylim=(0,1))
    # Plots.annotate!(p, [(0.02, 0.95, Plots.text(stats, 10, :left))]; subplot=2)

    return p
end
