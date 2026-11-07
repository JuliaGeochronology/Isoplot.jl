module PlotsExt

    using Isoplot
    using Plots: Shape, plot, plot!, plot3d!
    import Plots
    using Measurements
    import Isoplot: concordiacoords, concordialinecoeffs, concordiacurvedata, xtimerange, wmeanband, ellipsoid, calibrationlabels, calibrationlabel
    # export plot, plot!

    const PlotOrSubplot = Union{Plots.Plot, Plots.Subplot}
    Base.retry_load_extensions() 

    # Generic redirect from `plot` to `plot!` (only implement methods for `plot!` subsequently)
    Plots.plot(d::Union{Data,Vector{<:Data}}, args...; framestyle=:box, kwargs...) = plot!(plot(), d, args...; framestyle, kwargs...)
    Plots.plot(x, y::Union{Data,Vector{<:Data}}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)
    Plots.plot(x::Union{Data,Vector{<:Data}}, y::Union{Data,Vector{<:Data}}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)

    # Plot 2d uncertainty ellipses of any sort
    Plots.Shape(e::Ellipse{T}) where {T} = Shape{T,T}(e.x, e.y)
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), a::AbstractAnalysis, args...; kwargs...) = plot!(hdl, Ellipse(a), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), a::Vector{<:AbstractAnalysis}, args...; kwargs...) = plot!(hdl, Ellipse.(a), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), e::Ellipse, args...; kwargs...) = plot!(hdl, Shape(e), args...; kwargs...)
        @eval function Plots.plot!(hdl::($P), e::Vector{<:Ellipse}, args...; label="", color=1, kwargs...)
            for i in eachindex(e)
                plot!(hdl, Shape(e[i]), args...; label=(i==firstindex(e) ? label : ""), color, kwargs...)
            end
            return hdl
        end
    end

    # 3D analyses are plotted in 2D as the ellipse of any pair of variables `dims` (default 1 and 2),
    # or in 3D as an ellipsoid by `plot3d!`
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), a::Analysis3D, args...; dims=(1,2), kwargs...) = plot!(hdl, Ellipse(a, dims...), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), a::Vector{<:Analysis3D}, args...; dims=(1,2), kwargs...) = plot!(hdl, Ellipse.(a, dims...), args...; kwargs...)
        @eval Plots.plot3d!(hdl::($P), a::Analysis3D; color=1, alpha=0.3, label="", sigmalevel=2.7954, npoints=24, kwargs...) =
            plot!(hdl, ellipsoid(a; sigmalevel, npoints)...; seriestype=:surface, color=Plots.cgrad(fill(color isa Integer ? Plots.palette(:default)[color] : color, 2)), colorbar=false, alpha, label, kwargs...)
        @eval function Plots.plot3d!(hdl::($P), a::Vector{<:Analysis3D}; label="", color=1, kwargs...)
            for i in eachindex(a)
                plot3d!(hdl, a[i]; label=(i==firstindex(a) ? label : ""), color, kwargs...)
            end
            return hdl
        end
    end

    # Line fits are plotted in 2D by `plot!` (variable `dims` against the first, with 2σ uncertainty),
    # and in 3D by `plot3d!` (variables `dims` = (j, k) against the first, without uncertainty)
    for P in (Plots.Plot, Plots.Subplot)
        @eval function Plots.plot!(hdl::($P), lf::Union{YorkFit, McLeanFit}; npoints=25, dims::Integer=2, kwargs...)
            x = collect(range(Plots.xlims(hdl)..., length=npoints))
            y = Isoplot.line.(lf, x, dims)
            return plot!(hdl, x, value.(y); ribbon=2stdev.(y), kwargs...)
        end
        @eval function Plots.plot3d!(hdl::($P), lf::McLeanFit; npoints=25, dims=(2,3), kwargs...)
            x = collect(range(Plots.xlims(hdl)..., length=npoints))
            y, z = (value.(Isoplot.line.(lf, x, d)) for d in dims)
            return plot3d!(hdl, x, y, z; kwargs...)
        end
    end

    Plots.plot(calib::Isoplot.UPbSIMSCalibration; framestyle=:box, kwargs...) = plot!(plot(), calib; framestyle, kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval function Plots.plot!(hdl::($P), calib::Isoplot.UPbSIMSCalibration; 
                color = :auto,
                alpha = 0.75,
                xlabel=calibrationlabels(calib).xlabel,
                ylabel=calibrationlabels(calib).ylabel,
                kwargs...
            )
            plot!(hdl, calib.data; color=(color===:auto ? 1 : color), alpha, xlabel, ylabel, label="Data (95% CI)", kwargs...)
            plot!(hdl, calib.line; color=(color===:auto ? 2 : color), label=calibrationlabel(calib), kwargs...)
            return hdl
        end
    end

    ## --- Concordia lines and curves, in Wetherill space (`concordia*`) or Tera-Wasserburg space (`twconcordia*`)

    for (tw, line, curve) in ((false, :concordialine, :concordiacurve), (true, :twconcordialine, :twconcordiacurve))
        # Plot a line between two times, or between many pairs of times
        @eval function Isoplot.$line(t₀, t₁; framestyle=:box, kwargs...)
            xl = extrema(first(concordiacoords((first(t₀), first(t₁)), value(λ235U), value(λ238U); tw=$tw)))
            return Isoplot.$(Symbol(line, :!))(plot(xlims=xl), t₀, t₁; framestyle, kwargs...)
        end
        @eval Isoplot.$(Symbol(line, :!))(hdl::PlotOrSubplot, t₀, t₁; kwargs...) = plotconcordialine!(hdl, t₀, t₁, $tw; kwargs...)
        # Plot the concordia curve
        @eval Isoplot.$(Symbol(curve, :!))(hdl::PlotOrSubplot=Plots.current(); kwargs...) = plotconcordiacurve!(hdl, $tw; kwargs...)
    end

    # Points spanning the x limits of a plot, or just the part of them between `xmin` and `xmax` if `truncate`
    xrange(xl, xmin, xmax, truncate::Bool) = collect(range((truncate ? (max(first(xl), xmin), min(last(xl), xmax)) : xl)..., length=50))

    # Single line, with uncertainty from the decay constants
    function plotconcordialine!(hdl::PlotOrSubplot, t₀::Number, t₁::Number, tw::Bool; truncate::Bool=false, kwargs...)
        xl = Plots.xlims(hdl)
        intercept, slope, x₀, x₁ = concordialinecoeffs(t₀, t₁; tw)
        x = xrange(xl, min(value(x₀), value(x₁)), max(value(x₀), value(x₁)), truncate)
        y = intercept .+ slope .* x
        plot!(hdl, x, value.(y); ribbon=2stdev.(y), kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end
    # Many lines, plus their mean
    function plotconcordialine!(hdl::PlotOrSubplot, t₀::Collection, t₁::Collection, tw::Bool; truncate::Bool=false, label="", color=:black, alpha=0.05, kwargs...)
        xl = Plots.xlims(hdl)
        intercept, slope, x₀, x₁ = concordialinecoeffs(t₀, t₁, value(λ235U), value(λ238U); tw)
        x = xrange(xl, min(minimum(x₀), minimum(x₁)), max(maximum(x₀), maximum(x₁)), truncate)
        ys = [@. b + m * x for (m, b) in zip(slope, intercept)]
        plot!(hdl, x, ys; label="", color, alpha, kwargs...)
        plot!(hdl, x, sum(ys)./length(ys); label, color, alpha=1, kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end

    function plotconcordiacurve!(hdl::PlotOrSubplot, tw::Bool; timerange=:auto)
        xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
        tlim = timerange===:auto ? xtimerange(xl; tw) : timerange # Time range of current window
        c = concordiacurvedata(tlim; tw)
        Plots.plot!(hdl, Shape(c.envelope...), color=:black, alpha=0.15, label="") # Two-sigma concordia uncertainty
        Plots.plot!(hdl, c.curve..., color=:black, label="") # Concordia line
        Plots.plot!(hdl, c.schoene..., color=:black, style=:dash, label="") # Concordia line, using Schoene et al. decay constants
        Plots.xlims!(hdl, xl) # Ensure that figure size hasn't changed
        Plots.ylims!(hdl, yl)

        # Plot age markers with text labels
        xticks, yticks = c.ticks
        Plots.plot!(hdl, xticks, yticks, color=:black, seriestype=:scatter, ms=2, label="")
        xoffset = (xl[2]-xl[1])/200
        yoffset = (yl[2]-yl[1])/100
        visible = (xl[1]+8*xoffset) .< xticks .< xl[2]
        ticklabels = Plots.text.(string.(round.(c.tticks[visible], digits=1)), 10, :right)
        Plots.annotate!(hdl, xticks[visible].-xoffset, yticks[visible].+yoffset, ticklabels)
        return hdl
    end

    # Plot confidence intervals
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), y::Vector{<:CI}; mscolor=:auto, kwargs...) = plot!(hdl, y.|>c->c.mean; yerror=(y.|>c->c.mean-c.lower, y.|>c->c.upper-c.mean), mscolor, kwargs...)
        @eval Plots.plot!(hdl::($P), x, y::Vector{<:CI}; mscolor=:auto, kwargs...) = plot!(hdl, x, y.|>c->c.mean; yerror=(y.|>c->c.mean-c.lower, y.|>c->c.upper-c.mean), mscolor, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:CI}, y; mscolor=:auto, kwargs...) = plot!(hdl, x.|>c->c.mean, y; xerror=(x.|>c->c.mean-c.lower, x.|>c->c.upper-c.mean), mscolor, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:CI}, y::Vector{<:CI}; mscolor=:auto, kwargs...) = plot!(hdl, x.|>c->c.mean, y.|>c->c.mean; xerror=(x.|>c->c.mean-c.lower, x.|>c->c.upper-c.mean), yerror=(y.|>c->c.mean-c.lower, y.|>c->c.upper-c.mean), mscolor, kwargs...)
    end

    # Rank-order plots
    Isoplot.rankorder(args...; framestyle=:box, kwargs...) = rankorder!(plot(), args...; framestyle, kwargs...)
    Isoplot.rankorder!(h::PlotOrSubplot, data, sigma, i0::Number=0; kwargs...) = rankorder!(h, data .± sigma, i0; kwargs...)
    function Isoplot.rankorder!(h::PlotOrSubplot, data::Vector{<:Measurement}, i0::Number=0;
            scale=1,
            label="",
            mscolor=:auto,
            seriestype=:scatter,
            xticks=Float64[],
            wmean=false, # Overlay the weighted mean and its 2σ interval
            kwargs...
        )
        x = i0 .+ scale.*(1:length(data))
        ds = sort(data)
        plot!(h, x, value.(ds), yerror=2*stdev.(ds); label, mscolor, seriestype, xticks, kwargs...)
        if wmean
            m = wmeanband(data)
            plot!(h, [first(x)-scale/2, last(x)+scale/2], [m.μ, m.μ]; ribbon=m.σ, color=:black, fillalpha=0.15, label=m.label)
        end
        return h
    end
    function Isoplot.rankorder!(h::PlotOrSubplot, data::Vector{<:CI}, i0::Number=0;
            scale=1,
            label="",
            mscolor=:auto,
            seriestype=:scatter,
            xticks=Float64[],
            wmean=false,
            kwargs...
        )
        wmean && throw(ArgumentError("`wmean` is only supported for `Measurement` data"))
        x = i0 .+ scale.*(1:length(data))
        y = sort(data)
        plot!(h, x, y.|>c->c.mean; yerror=(y.|>c->c.mean-c.lower, y.|>c->c.upper-c.mean), label, mscolor, seriestype, xticks, kwargs...)
    end

end