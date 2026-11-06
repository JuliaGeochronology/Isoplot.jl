module PlotsExt

    using Isoplot
    using Plots: Shape, plot, plot!
    import Plots
    using Measurements
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

    for P in (Plots.Plot, Plots.Subplot)
        @eval function Plots.plot!(hdl::($P), yf::YorkFit; npoints=25, kwargs...)
            x = range(Plots.xlims(hdl)..., length=npoints)
            y = Isoplot.line.(yf, x)
            plot!(hdl, x, value.(y); ribbon=2stdev.(y), kwargs...)
            return hdl
        end
    end

    Plots.plot(calib::Isoplot.UPbSIMSCalibration; framestyle=:box, kwargs...) = plot!(plot(), calib; framestyle, kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval function Plots.plot!(hdl::($P), calib::Isoplot.UPbSIMSCalibration; 
                color = :auto,
                alpha = 0.75,
                xlabel="$(Isoplot.prettify(:Pb206)) / $(Isoplot.prettify(calib.U)) RSF",
                ylabel="$(Isoplot.prettify(calib.numerator)) / $(Isoplot.prettify(calib.denominator))",
                kwargs...
            )
            plot!(hdl, calib.data; color=(color===:auto ? 1 : color), alpha, xlabel, ylabel, label="Data (95% CI)", kwargs...)
            plot!(hdl, calib.line; color=(color===:auto ? 2 : color), label="York fit (N=$(length(calib.data)), MSWD=$(round(calib.line.mswd, sigdigits=3)))", kwargs...)
            return hdl
        end
    end

    # Plot a line between two times in Wetherill Concordia space
    Isoplot.concordialine(t₀, t₁; framestyle=:box, kwargs...) = concordialine!(plot(xlims=ratio.((first(t₀), first(t₁)), value(λ235U))), t₀, t₁; framestyle, kwargs...)
    function Isoplot.concordialine!(hdl::PlotOrSubplot, t₀::Number, t₁::Number; truncate::Bool=false, kwargs...)
        xl = Plots.xlims(hdl)
        r75₀ = ratio(t₀, value(λ235U))
        r68₀ = ratio(t₀, value(λ238U))
        r75₁ = ratio(t₁, value(λ235U))
        r68₁ = ratio(t₁, value(λ238U))
        slope = (r68₁-r68₀)/(r75₁-r75₀)
        intercept = r68₀ - r75₀*slope
        x = if truncate
            xmin = max(first(xl), value(r75₀))
            xmax = min(last(xl), value(r75₁))
            range(xmin, xmax, length=50)
        else
            range(xl..., length=50)
        end
        y = intercept .+ slope .* x
        plot!(hdl, x, value.(y); ribbon=2stdev.(y), kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end
    function Isoplot.concordialine!(hdl::PlotOrSubplot, t₀::Collection, t₁::Collection; truncate::Bool=false, label="", color=:black, alpha=0.05, kwargs...)
        xl = Plots.xlims(hdl)
        r75₀ = ratio.(t₀, value(λ235U))
        r68₀ = ratio.(t₀, value(λ238U))
        r75₁ = ratio.(t₁, value(λ235U))
        r68₁ = ratio.(t₁, value(λ238U))
        slope = @. (r68₁-r68₀)/(r75₁-r75₀)
        intercept = @. r68₀ - r75₀*slope
        x = if truncate
            xmin = max(first(xl), minimum(r75₀))
            xmax = min(last(xl), maximum(r75₁))
            collect(range(xmin, xmax, length=50))
        else
            collect(range(xl..., length=50))
        end
        y(slope, intercept) = @. intercept + slope * x
        ys = y.(slope, intercept)
        plot!(hdl, x, ys; label="", color, alpha, kwargs...)
        plot!(hdl, x, sum(ys)./length(ys); label, color, alpha=1, kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end

    # Plot the Wetherill Concordia curve
    function Isoplot.concordiacurve!(hdl::PlotOrSubplot=Plots.current(); timerange=:auto)
        # Uncertainty of 235 decay constant relative to the 238 decay constant
        σₜ = value(λ235U_jaffey) .* sqrt((stdev(λ238U)/value(λ238U)).^2 + (stdev(λ235U_jaffey)/value(λ235U_jaffey)).^2) # 1/Years

        # Plot the concordia curve
        xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
        tlim = if timerange===:auto
             age.(max.(xl, 0.0), value(λ235U_jaffey)) # Calculate time range of current window
        else
            timerange
        end
        dt = tlim[2] - tlim[1]
        tmin = max(tlim[1]-0.1dt, 0.0)
        tmax = tlim[2]+0.1dt
        t = range(tmin, tmax, length=1000) # Time vector, including padding
        r75t = ratio.(t, value(λ235U_jaffey)) # X axis values
        r68t = ratio.(t, value(λ238U)) # Y axis values
        x = [ratio.(t, value(λ235U_jaffey)-σₜ*2); reverse(ratio.(t, value(λ235U_jaffey)+σₜ*2))]
        y = [r68t; reverse(r68t)]
        Plots.plot!(hdl, Shape(x,y), color=:black, alpha=0.15, label="") # Two-sigma concordia uncertainty
        Plots.plot!(hdl, r75t, r68t, color=:black, label="") # Concordia line

        r75t_Schoene = ratio.(t, value(λ235U)) # X axis values
        Plots.plot!(hdl, r75t_Schoene,r68t,color=:black,style=:dash, label="") # Concordia line

        Plots.xlims!(hdl, xl) # Ensure that figure size hasn't changed
        Plots.ylims!(hdl, yl)

        # Calculate desired range of age markers
        scale = floor(log10(tlim[2]-tlim[1])) # Characteristic timescale (order of magnitude)
        trange = round.(tlim./10.0^scale) # Minimum and maximum time to a round number
        majorstep = 0.5
        tticks = (trange[1]:majorstep:trange[2]).*10.0^scale # Time ticks, to a round number
        r75tticks = ratio.(tticks, value(λ235U_jaffey)) # X axis values
        r68tticks = ratio.(tticks, value(λ238U)) # Y axis values

        # Plot age markers with text labels
        Plots.plot!(hdl, r75tticks,r68tticks, color=:black, seriestype=:scatter, ms=2, label="")
        xoffset = (xl[2]-xl[1])/200
        yoffset = (yl[2]-yl[1])/100
        t = (xl[1]+8*xoffset) .< r75tticks .< xl[2]
        ticklabels = Plots.text.(string.(round.(tticks[t], digits=1)), 10, :right)
        Plots.annotate!(hdl, r75tticks[t].-xoffset,r68tticks[t].+yoffset,ticklabels)

        return hdl
    end

    ## --- Tera-Wasserburg concordia (²³⁸U/²⁰⁶Pb vs. ²⁰⁷Pb/²⁰⁶Pb)

    # Tera-Wasserburg coordinates for a given time, from the Wetherill ratios
    # (Works with Numbers or Measurements)
    function twcoords(t, λ235, λ238)
        r75 = ratio(t, λ235)
        r68 = ratio(t, λ238)
        return 1/r68, r75/r68/value(Isoplot.U238_U235)
    end

    # Plot a line between two times in Tera-Wasserburg Concordia space
    # (A straight line in Wetherill space is also a straight line in Tera-Wasserburg space)
    function Isoplot.twconcordialine(t₀, t₁; framestyle=:box, kwargs...)
        xl = extrema(first(twcoords(t, value(λ235U), value(λ238U))) for t in (first(t₀), first(t₁)))
        return twconcordialine!(plot(xlims=xl), t₀, t₁; framestyle, kwargs...)
    end
    function Isoplot.twconcordialine!(hdl::PlotOrSubplot, t₀::Number, t₁::Number; truncate::Bool=false, kwargs...)
        xl = Plots.xlims(hdl)
        x₀, y₀ = twcoords(t₀, λ235U, λ238U)
        x₁, y₁ = twcoords(t₁, λ235U, λ238U)
        slope = (y₁-y₀)/(x₁-x₀)
        intercept = y₀ - x₀*slope
        x = if truncate
            xmin = max(first(xl), min(value(x₀), value(x₁)))
            xmax = min(last(xl), max(value(x₀), value(x₁)))
            collect(range(xmin, xmax, length=50))
        else
            collect(range(xl..., length=50))
        end
        y = intercept .+ slope .* x
        plot!(hdl, x, value.(y); ribbon=2stdev.(y), kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end
    function Isoplot.twconcordialine!(hdl::PlotOrSubplot, t₀::Collection, t₁::Collection; truncate::Bool=false, label="", color=:black, alpha=0.05, kwargs...)
        xl = Plots.xlims(hdl)
        x₀, y₀ = twcoords.(t₀, value(λ235U), value(λ238U)) |> v -> (first.(v), last.(v))
        x₁, y₁ = twcoords.(t₁, value(λ235U), value(λ238U)) |> v -> (first.(v), last.(v))
        slope = @. (y₁-y₀)/(x₁-x₀)
        intercept = @. y₀ - x₀*slope
        x = if truncate
            xmin = max(first(xl), min(minimum(x₀), minimum(x₁)))
            xmax = min(last(xl), max(maximum(x₀), maximum(x₁)))
            collect(range(xmin, xmax, length=50))
        else
            collect(range(xl..., length=50))
        end
        ys = [@. b + m * x for (m, b) in zip(slope, intercept)]
        plot!(hdl, x, ys; label="", color, alpha, kwargs...)
        plot!(hdl, x, sum(ys)./length(ys); label, color, alpha=1, kwargs...)
        Plots.xlims!(hdl, xl)
        return hdl
    end

    # Plot the Tera-Wasserburg Concordia curve
    function Isoplot.twconcordiacurve!(hdl::PlotOrSubplot=Plots.current(); timerange=:auto)
        # Uncertainty of 235 decay constant relative to the 238 decay constant
        σₜ = value(λ235U_jaffey) .* sqrt((stdev(λ238U)/value(λ238U)).^2 + (stdev(λ235U_jaffey)/value(λ235U_jaffey)).^2) # 1/Years
        λ235, λ238 = value(λ235U_jaffey), value(λ238U)

        xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
        tlim = if timerange===:auto
            tmax = 4567.0 # Don't extrapolate beyond the age of the Earth (Ma)
            xmin = max(first(xl), 1/ratio(tmax, λ238)) # Time decreases with x
            xmax = max(last(xl), xmin)
            (age(1/xmax, λ238), age(1/xmin, λ238)) # Calculate time range of current window
        else
            timerange
        end
        dt = tlim[2] - tlim[1]
        tmin = max(tlim[1]-0.1dt, 1e-3dt) # Tera-Wasserburg space diverges at t=0
        tmax = tlim[2]+0.1dt
        t = range(tmin, tmax, length=1000) # Time vector, including padding
        x, y = twcoords.(t, λ235, λ238) |> v -> (first.(v), last.(v))
        _, ylo = twcoords.(t, λ235-2σₜ, λ238) |> v -> (first.(v), last.(v))
        _, yhi = twcoords.(t, λ235+2σₜ, λ238) |> v -> (first.(v), last.(v))
        Plots.plot!(hdl, Shape([x; reverse(x)], [ylo; reverse(yhi)]), color=:black, alpha=0.15, label="") # Two-sigma concordia uncertainty
        Plots.plot!(hdl, x, y, color=:black, label="") # Concordia line

        xs, ys = twcoords.(t, value(λ235U), λ238) |> v -> (first.(v), last.(v))
        Plots.plot!(hdl, xs, ys, color=:black, style=:dash, label="") # Concordia line, using Schoene et al. decay constants

        Plots.xlims!(hdl, xl) # Ensure that figure size hasn't changed
        Plots.ylims!(hdl, yl)

        # Calculate desired range of age markers
        scale = floor(log10(tlim[2]-tlim[1])) # Characteristic timescale (order of magnitude)
        trange = round.(tlim./10.0^scale) # Minimum and maximum time to a round number
        majorstep = 0.5
        tticks = collect((trange[1]:majorstep:trange[2]).*10.0^scale) # Time ticks, to a round number
        filter!(τ -> tmin <= τ <= tmax, tticks)
        xticks, yticks = twcoords.(tticks, λ235, λ238) |> v -> (first.(v), last.(v))

        # Plot age markers with text labels
        Plots.plot!(hdl, xticks, yticks, color=:black, seriestype=:scatter, ms=2, label="")
        xoffset = (xl[2]-xl[1])/200
        yoffset = (yl[2]-yl[1])/100
        visible = (xl[1]+8*xoffset) .< xticks .< xl[2]
        ticklabels = Plots.text.(string.(round.(tticks[visible], digits=1)), 10, :right)
        Plots.annotate!(hdl, xticks[visible].-xoffset, yticks[visible].+yoffset, ticklabels)

        return hdl
    end

    # Plot confidence intervals
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), y::Vector{<:CI}; mscolor=:auto, kwargs...) = plot!(hdl, y.|>c->c.mean; yerror=(y.|>c->c.upper-c.mean, y.|>c->c.mean-c.lower), mscolor, kwargs...)
        @eval Plots.plot!(hdl::($P), x, y::Vector{<:CI}; mscolor=:auto, kwargs...) = plot!(hdl, x, y.|>c->c.mean; yerror=(y.|>c->c.upper-c.mean, y.|>c->c.mean-c.lower), mscolor, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:CI}, y; mscolor=:auto, kwargs...) = plot!(hdl, x.|>c->c.mean, y; xerror=(x.|>c->c.upper-c.mean, x.|>c->c.mean-c.lower), mscolor, kwargs...)
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
            kwargs...
        )
        x = i0 .+ scale.*(1:length(data))
        ds = sort(data)
        plot!(h, x, value.(ds), yerror=2*stdev.(ds); label, mscolor, seriestype, xticks, kwargs...)
    end
    function Isoplot.rankorder!(h::PlotOrSubplot, data::Vector{<:CI}, i0::Number=0;
            scale=1,
            label="",
            mscolor=:auto,
            seriestype=:scatter,
            xticks=Float64[],
            kwargs...
        )
        x = i0 .+ scale.*(1:length(data))
        y = sort(data)
        plot!(h, x, y.|>c->c.mean; yerror=(y.|>c->c.upper-c.mean, y.|>c->c.mean-c.lower), label, mscolor, seriestype, xticks, kwargs...)
    end

end