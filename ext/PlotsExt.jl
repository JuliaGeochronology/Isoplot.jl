module PlotsExt

    using Isoplot
    using Plots: Shape, plot, plot!
    import Plots
    using Measurements
    # export plot, plot!

    const PlotOrSubplot = Union{Plots.Plot, Plots.Subplot}
    Base.retry_load_extensions() 
    # Plot 2d uncertainty ellipses of any sort
    Plots.Shape(e::Ellipse{T}) where {T} = Shape{T,T}(e.x, e.y)

    Plots.plot(e::Union{Data,Vector{<:Data}}, args...; kwargs...) = plot!(plot(), e, args...; kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), a::Analysis, args...; kwargs...) = plot!(hdl, Ellipse(a), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), a::Vector{<:Analysis}, args...; kwargs...) = plot!(hdl, Ellipse.(a), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), e::Ellipse, args...; kwargs...) = plot!(hdl, Shape(e), args...; kwargs...)
        @eval Plots.plot!(hdl::($P), e::Vector{<:Ellipse}, args...; kwargs...) = plot!(hdl, Shape.(e), args...; kwargs...)
    end

    # Plot a line between two times in Wetherill Concordia space
    Isoplot.concordialine(t₀, t₁; framestyle=:box, kwargs...) = concordialine!(plot(xlims=ratio.((first(t₀), first(t₁)), λ235U.val)), t₀, t₁; framestyle, kwargs...)
    function Isoplot.concordialine!(hdl::PlotOrSubplot, t₀::Number, t₁::Number; truncate::Bool=false, kwargs...)
        xl = Plots.xlims(hdl)
        r75₀ = ratio(t₀, λ235U.val)
        r68₀ = ratio(t₀, λ238U.val)
        r75₁ = ratio(t₁, λ235U.val)
        r68₁ = ratio(t₁, λ238U.val)
        slope = (r68₁-r68₀)/(r75₁-r75₀)
        intercept = r68₀ - r75₀*slope
        x = if truncate
            xmin = max(first(xl), val(r75₀))
            xmax = min(last(xl), val(r75₁))
            range(xmin, xmax, length=50)
        else
            range(xl..., length=50)
        end
        y = intercept .+ slope .* x
        plot!(hdl, x, val.(y); ribbon=err.(y), kwargs...)
        Plots.xlims!(hdl, xl)
    end
    function Isoplot.concordialine!(hdl::PlotOrSubplot, t₀::Collection, t₁::Collection; truncate::Bool=false, label="", color=:black, alpha=0.05, kwargs...)
        xl = Plots.xlims(hdl)
        r75₀ = ratio.(t₀, λ235U.val)
        r68₀ = ratio.(t₀, λ238U.val)
        r75₁ = ratio.(t₁, λ235U.val)
        r68₁ = ratio.(t₁, λ238U.val)
        slope = @. (r68₁-r68₀)/(r75₁-r75₀)
        intercept = @. r68₀ - r75₀*slope
        x = if truncate
            xmin = max(first(xl), vminimum(r75₀))
            xmax = min(last(xl), vmaximum(r75₁))
            collect(range(xmin, xmax, length=50))
        else
            collect(range(xl..., length=50))
        end
        y(slope, intercept) = @. intercept + slope * x
        ys = y.(slope, intercept)
        plot!(hdl, x, ys; label="", color, alpha, kwargs...)
        plot!(hdl, x, sum(ys)./length(ys); label, color, alpha=1, kwargs...)
        Plots.xlims!(hdl, xl)
    end

    # Plot the Wetherill Concordia curve
    function Isoplot.concordiacurve!(hdl::PlotOrSubplot=Plots.current())
        # Uncertainty of 235 decay constant relative to the 238 decay constant
        σₜ = λ235U_jaffey.val .* sqrt((λ238U.err/λ238U.val).^2 + (λ235U_jaffey.err/λ235U_jaffey.val).^2) # 1/Years

        # Plot the concordia curve
        xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
        tlim = age.(max.(xl, 0.0), λ235U_jaffey.val) # Calculate time range of current window
        dt = tlim[2] - tlim[1]
        tmin = max(tlim[1]-0.1dt, 0.0)
        tmax = tlim[2]+0.1dt
        t = range(tmin, tmax, length=1000) # Time vector, including padding
        r75t = ratio.(t, λ235U_jaffey.val) # X axis values
        r68t = ratio.(t, λ238U.val) # Y axis values
        x = [ratio.(t, λ235U_jaffey.val-σₜ*2); reverse(ratio.(t, λ235U_jaffey.val+σₜ*2))]
        y = [r68t; reverse(r68t)]
        Plots.plot!(hdl, Shape(x,y), color=:black, alpha=0.15, label="") # Two-sigma concordia uncertainty
        Plots.plot!(hdl, r75t, r68t, color=:black, label="") # Concordia line

        r75t_Schoene = ratio.(t, λ235U.val) # X axis values
        Plots.plot!(hdl, r75t_Schoene,r68t,color=:black,style=:dash, label="") # Concordia line

        Plots.xlims!(hdl, xl) # Ensure that figure size hasn't changed
        Plots.ylims!(hdl, yl)

        # Calculate desired range of age markers
        scale = floor(log10(tlim[2]-tlim[1])) # Characteristic timescale (order of magnitude)
        trange = round.(tlim./10.0^scale) # Minimum and maximum time to a round number
        majorstep = 0.5
        tticks = (trange[1]:majorstep:trange[2]).*10.0^scale # Time ticks, to a round number
        r75tticks = ratio.(tticks, λ235U_jaffey.val) # X axis values
        r68tticks = ratio.(tticks, λ238U.val) # Y axis values

        # Plot age markers with text labels
        Plots.plot!(hdl, r75tticks,r68tticks, color=:black, seriestype=:scatter, ms=2, label="")
        xoffset = (xl[2]-xl[1])/200
        yoffset = (yl[2]-yl[1])/100
        t = (xl[1]+8*xoffset) .< r75tticks .< xl[2]
        ticklabels = Plots.text.(string.(round.(Int, tticks[t])), 10, :right)
        Plots.annotate!(hdl, r75tticks[t].-xoffset,r68tticks[t].+yoffset,ticklabels)

        return hdl
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
        plot!(h, x, sort(data); label, mscolor, seriestype, xticks, kwargs...)
    end
end