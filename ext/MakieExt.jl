module MakieExt

    using Isoplot
    using Makie
    using Measurements
    import Isoplot: concordialine, concordialine!, concordiacurve, concordiacurve!

    ## --- Ellipses and analyses

    Makie.convert_arguments(P::Type{<:Poly}, e::Ellipse) = convert_arguments(P, as_points(e.x, e.y))
    Makie.convert_arguments(P::Type{<:Poly}, a::AbstractAnalysis) = convert_arguments(P, Ellipse(a))
    # Many ellipses at once, as a vector of polygons
    Makie.convert_arguments(::Type{<:Poly}, es::AbstractVector{<:Ellipse}) = (Makie.GeometryBasics.Polygon.(as_points.(getfield.(es, :x), getfield.(es, :y))),)
    Makie.convert_arguments(P::Type{<:Poly}, as::AbstractVector{<:AbstractAnalysis}) = convert_arguments(P, Ellipse.(as))
    Makie.plottype(::Ellipse) = Poly
    Makie.plottype(::AbstractAnalysis) = Poly
    Makie.plottype(::AbstractVector{<:Ellipse}) = Poly
    Makie.plottype(::AbstractVector{<:AbstractAnalysis}) = Poly

    function as_points(v1, v2)
        length(v1) == length(v2) || throw(DimensionMismatch("Arguments must be the same length"))
        return Point2f[Point2f(v1[i], v2[i]) for i in 1:length(v1)] # (Not `zip`, which would make a statically-sized result for SVectors)
    end

    ## --- Helpers

    # Current limits of the active axis, as ((xmin, xmax), (ymin, ymax))
    function axislimits()
        ax = Makie.current_axis()
        ax isa Makie.Axis || throw(ArgumentError("No current Axis to take limits from; pass the time range explicitly"))
        Makie.reset_limits!(ax)
        r = ax.finallimits[]
        return (r.origin[1], r.origin[1] + r.widths[1]), (r.origin[2], r.origin[2] + r.widths[2])
    end

    # Convert a time and decay constants to coordinates in Wetherill (²⁰⁷Pb/²³⁵U, ²⁰⁶Pb/²³⁸U)
    # or Tera-Wasserburg (²³⁸U/²⁰⁶Pb, ²⁰⁷Pb/²⁰⁶Pb) space. Works on Numbers or Measurements
    function concordiacoords(t, λ235, λ238, concordiatype::Symbol)
        r75 = ratio(t, λ235)
        r68 = ratio(t, λ238)
        if concordiatype === :Wetherill
            return r75, r68
        elseif concordiatype === :TeraWasserburg
            return 1/r68, r75/r68/value(Isoplot.U238_U235)
        else
            throw(ArgumentError("concordiaType must be :Wetherill or :TeraWasserburg"))
        end
    end

    # Vectorized version, returning a vector of x values and a vector of y values
    function concordiacoords(t::AbstractArray, λ235, λ238, concordiatype::Symbol)
        xy = concordiacoords.(t, λ235, λ238, concordiatype)
        return first.(xy), last.(xy)
    end

    # Convert a position along the x axis back into a time
    function timefromx(x, concordiatype::Symbol)
        if concordiatype === :Wetherill
            return age(x, value(λ235U_jaffey))
        elseif concordiatype === :TeraWasserburg
            return age(1/x, value(λ238U))
        else
            throw(ArgumentError("concordiaType must be :Wetherill or :TeraWasserburg"))
        end
    end

    # Time range corresponding to the x limits of the current axis
    function autotimerange(concordiatype)
        xl, _ = axislimits()
        tmax = 4567.0 # Don't extrapolate beyond the age of the Earth (Ma)
        tl = if concordiatype === :Wetherill
            timefromx.(max.(xl, 0.0), :Wetherill)
        else # Time decreases with x
            sort(timefromx.(max.(xl, 1/ratio(tmax, value(λ238U))), :TeraWasserburg))
        end
        return clamp(tl[1], 0.0, tmax), clamp(tl[2], 0.0, tmax)
    end

    ## --- Lines between two points on the concordia

    # Plot a line between two times in Concordia space
    @recipe(ConcordiaLine, t₀, t₁) do scene
        Attributes(
            concordiaType = :Wetherill,
            color = :black,
            linewidth = 1,
            linestyle = :solid,
            errorcolor = (:black, 0.15),
            alpha = 0.05, # Transparency of individual lines, when plotting many at once
            label = nothing,
        )
    end

    function Makie.plot!(cline::ConcordiaLine)
        type = cline.concordiaType[]
        t₀, t₁ = cline.t₀[], cline.t₁[]
        x₀, y₀ = concordiacoords(t₀, λ235U, λ238U, type)
        x₁, y₁ = concordiacoords(t₁, λ235U, λ238U, type)
        slope = @. (y₁ - y₀)/(x₁ - x₀)
        intercept = @. y₀ - x₀*slope
        if t₀ isa Number && t₁ isa Number
            # Single line, with decay constant uncertainty
            x = collect(range(value(x₀), value(x₁), length=50))
            y = intercept .+ slope .* x
            ŷ, σ = value.(y), 2stdev.(y)
            band = [as_points(x, ŷ .+ σ); as_points(reverse(x), reverse(ŷ .- σ))]
            poly!(cline, band; color=cline.errorcolor, strokewidth=0)
            lines!(cline, x, ŷ; color=cline.color, linewidth=cline.linewidth, linestyle=cline.linestyle, label=cline.label)
        else
            # Many lines, plus their mean
            xs = collect(range(minimum(value, x₀), maximum(value, x₁), length=50))
            ys = [value.(b .+ m .* xs) for (m, b) in zip(slope, intercept)]
            for (i, y) in enumerate(ys)
                lines!(cline, xs, y; color=cline.color, alpha=cline.alpha, linewidth=cline.linewidth, linestyle=cline.linestyle)
            end
            lines!(cline, xs, sum(ys)./length(ys); color=cline.color, linewidth=cline.linewidth, linestyle=cline.linestyle, label=cline.label)
        end
        return cline
    end

    # Show only a simple line in legends, rather than all the child plots
    Makie.legendelements(p::ConcordiaLine, legend) = [LineElement(color=p.color[], linewidth=p.linewidth[], linestyle=p.linestyle[])]

    ## --- The concordia curve itself

    @recipe(ConcordiaCurve, t₀, t₁) do scene
        Attributes(
            concordiaType = :Wetherill,
            color = :black,
            linewidth = 1,
            errorcolor = (:black, 0.15),
            agefontsize = 10,
            ageticksize = 5,
            agetickcolor = :black,
            label = nothing,
        )
    end

    Makie.legendelements(p::ConcordiaCurve, legend) = [LineElement(color=p.color[], linewidth=p.linewidth[])]

    function Makie.plot!(ccurve::ConcordiaCurve)
        type = ccurve.concordiaType[]
        type in (:Wetherill, :TeraWasserburg) || throw(ArgumentError("concordiaType must be :Wetherill or :TeraWasserburg"))

        # Uncertainty of 235 decay constant relative to the 238 decay constant
        σₜ = value(λ235U_jaffey) .* sqrt((stdev(λ238U)/value(λ238U)).^2 + (stdev(λ235U_jaffey)/value(λ235U_jaffey)).^2) # 1/Years

        # Time range, with padding
        tlim = (ccurve.t₀[], ccurve.t₁[])
        dt = tlim[2] - tlim[1]
        tmin = max(tlim[1] - 0.1dt, 0.0)
        tmax = tlim[2] + 0.1dt
        t = range(tmin, tmax, length=1000)
        # Tera-Wasserburg space diverges at t=0, so start slightly after
        type === :TeraWasserburg && (t = range(max(tmin, 1e-3dt), tmax, length=1000))

        λ235 = value(λ235U_jaffey)
        λ238 = value(λ238U)
        x, y = concordiacoords(t, λ235, λ238, type)
        # Two-sigma uncertainty in the 235 decay constant
        xlo, ylo = concordiacoords(t, λ235 - 2σₜ, λ238, type)
        xhi, yhi = concordiacoords(t, λ235 + 2σₜ, λ238, type)
        if type === :Wetherill # Uncertainty is in x
            band = as_points([xlo; reverse(xhi)], [y; reverse(y)])
        else # Uncertainty is in y
            band = as_points([x; reverse(x)], [ylo; reverse(yhi)])
        end
        poly!(ccurve, band; color=ccurve.errorcolor, strokewidth=0)
        lines!(ccurve, x, y; color=ccurve.color, linewidth=ccurve.linewidth, label=ccurve.label)

        # Same curve using the Schoene et al. decay constants, dashed
        xs, ys = concordiacoords(t, value(λ235U), λ238, type)
        lines!(ccurve, xs, ys; color=ccurve.color, linewidth=ccurve.linewidth, linestyle=:dash)

        # Age markers at round numbers
        scale = floor(log10(tlim[2] - tlim[1])) # Characteristic timescale (order of magnitude)
        trange = round.(tlim ./ 10.0^scale) # Minimum and maximum time to a round number
        majorstep = 0.5
        tticks = collect((trange[1]:majorstep:trange[2]) .* 10.0^scale)
        filter!(τ -> max(tmin, 1e-3dt) <= τ <= tmax, tticks) # Only within the plotted curve (and Tera-Wasserburg is undefined at t=0)
        tickx, ticky = concordiacoords(tticks, λ235, λ238, type)
        scatter!(ccurve, tickx, ticky; markersize=ccurve.ageticksize, color=ccurve.agetickcolor)

        # Label the ticks that fall within the visible part of the curve
        xoffset = (maximum(x) - minimum(x))/200
        yoffset = (maximum(y) - minimum(y))/100
        visible = (minimum(x) + 8*xoffset) .< tickx .< maximum(x) # Labels sit to the left of their ticks
        ticklabels = string.(round.(tticks[visible], digits=1))
        text!(ccurve, tickx[visible] .- xoffset, ticky[visible] .+ yoffset; text=ticklabels, fontsize=ccurve.agefontsize, align=(:right, :bottom))

        return ccurve
    end

    # With no explicit time range, use the limits of the given axis
    concordiacurve!(ax::Makie.Axis; concordiaType=:Wetherill, timerange=nothing, kwargs...) = begin
        t₀, t₁ = timerange === nothing ? autotimerange(concordiaType) : timerange
        concordiacurve!(ax, t₀, t₁; concordiaType, kwargs...)
    end

    ## --- York fits

    @recipe(YorkLine, yf) do scene
        Attributes(
            color = Makie.wong_colors()[2],
            linewidth = 2,
            bandalpha = 0.3,
            npoints = 25,
            xlims = nothing, # Defaults to the limits of the current axis
            label = nothing,
        )
    end
    Makie.convert_arguments(::Type{<:YorkLine}, yf::YorkFit) = (yf,)
    Makie.plottype(::YorkFit) = YorkLine

    function Makie.plot!(p::YorkLine)
        yf = p.yf[]
        xl = p.xlims[] === nothing ? first(axislimits()) : p.xlims[]
        x = range(xl..., length=p.npoints[])
        y = Isoplot.line.(yf, x)
        ŷ, σ = value.(y), 2stdev.(y)
        band!(p, x, ŷ .- σ, ŷ .+ σ; color=(p.color[], p.bandalpha[]))
        lines!(p, x, ŷ; color=p.color, linewidth=p.linewidth, label=p.label)
        return p
    end

    ## --- Rank-order plots

    @recipe(RankOrderPlot, data) do scene
        Attributes(
            color = Makie.wong_colors()[1],
            markersize = 9,
            whiskerwidth = 0,
            scale = 1,
            offset = 0,
            label = nothing,
        )
    end
    Makie.convert_arguments(::Type{<:RankOrderPlot}, data::AbstractVector{<:Measurement}) = (data,)
    Makie.convert_arguments(::Type{<:RankOrderPlot}, data::AbstractVector{<:CI}) = (data,)
    Makie.convert_arguments(::Type{<:RankOrderPlot}, data::AbstractVector, sigma) = (data .± sigma,)

    # Methods are restricted to Makie axes so as not to collide with the Plots.jl extension
    Isoplot.rankorder!(ax::Makie.AbstractAxis, args...; kwargs...) = rankorderplot!(ax, args...; kwargs...)

    function Makie.plot!(p::RankOrderPlot)
        data = p.data[]
        n = length(data)
        x = p.offset[] .+ p.scale[] .* (1:n)
        if data isa AbstractVector{<:Measurement}
            ds = sort(data)
            μ, lo, hi = value.(ds), 2stdev.(ds), 2stdev.(ds)
        else
            ds = sort(data)
            μ = [c.mean for c in ds]
            lo = [c.mean - c.lower for c in ds]
            hi = [c.upper - c.mean for c in ds]
        end
        errorbars!(p, x, μ, lo, hi; color=p.color, whiskerwidth=p.whiskerwidth)
        scatter!(p, x, μ; color=p.color, markersize=p.markersize, label=p.label)
        return p
    end

end
