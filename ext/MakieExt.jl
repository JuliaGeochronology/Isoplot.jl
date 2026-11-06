module MakieExt

    using Isoplot
    using Makie
    using Measurements
    import Isoplot: concordialine, concordialine!, concordiacurve, concordiacurve!
    import Isoplot: concordiacoords, concordialinecoeffs, concordiacurvedata, xtimerange

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

    # Whether the `concordiaType` is Tera-Wasserburg (otherwise Wetherill)
    function istw(concordiatype::Symbol)
        concordiatype in (:Wetherill, :TeraWasserburg) || throw(ArgumentError("concordiaType must be :Wetherill or :TeraWasserburg"))
        return concordiatype === :TeraWasserburg
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
        tw = istw(cline.concordiaType[])
        intercept, slope, x₀, x₁ = concordialinecoeffs(cline.t₀[], cline.t₁[]; tw)
        if cline.t₀[] isa Number && cline.t₁[] isa Number
            # Single line, with decay constant uncertainty
            x = collect(range(value(x₀), value(x₁), length=50))
            y = intercept .+ slope .* x
            ŷ, σ = value.(y), 2stdev.(y)
            poly!(cline, [as_points(x, ŷ .+ σ); as_points(reverse(x), reverse(ŷ .- σ))]; color=cline.errorcolor, strokewidth=0)
            lines!(cline, x, ŷ; color=cline.color, linewidth=cline.linewidth, linestyle=cline.linestyle, label=cline.label)
        else
            # Many lines, plus their mean
            x = collect(range(minimum(value, x₀), maximum(value, x₁), length=50))
            ys = [value.(b .+ m .* x) for (m, b) in zip(slope, intercept)]
            for y in ys
                lines!(cline, x, y; color=cline.color, alpha=cline.alpha, linewidth=cline.linewidth, linestyle=cline.linestyle)
            end
            lines!(cline, x, sum(ys)./length(ys); color=cline.color, linewidth=cline.linewidth, linestyle=cline.linestyle, label=cline.label)
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
        tw = istw(ccurve.concordiaType[])
        tlim = (ccurve.t₀[], ccurve.t₁[])
        c = concordiacurvedata(tlim; tw)
        poly!(ccurve, as_points(c.envelope...); color=ccurve.errorcolor, strokewidth=0) # Two-sigma concordia uncertainty
        lines!(ccurve, c.curve...; color=ccurve.color, linewidth=ccurve.linewidth, label=ccurve.label)
        lines!(ccurve, c.schoene...; color=ccurve.color, linewidth=ccurve.linewidth, linestyle=:dash) # Using Schoene et al. decay constants

        # Age markers at round numbers, labeled where they fall within the visible part of the curve
        xticks, yticks = c.ticks
        scatter!(ccurve, xticks, yticks; markersize=ccurve.ageticksize, color=ccurve.agetickcolor)
        x, y = c.curve
        xoffset = (maximum(x) - minimum(x))/200
        yoffset = (maximum(y) - minimum(y))/100
        visible = (minimum(x) + 8*xoffset) .< xticks .< maximum(x) # Labels sit to the left of their ticks
        ticklabels = string.(round.(c.tticks[visible], digits=1))
        text!(ccurve, xticks[visible] .- xoffset, yticks[visible] .+ yoffset; text=ticklabels, fontsize=ccurve.agefontsize, align=(:right, :bottom))
        return ccurve
    end

    # With no explicit time range, use the limits of the given axis
    concordiacurve!(ax::Makie.Axis; concordiaType=:Wetherill, timerange=nothing, kwargs...) = begin
        t₀, t₁ = timerange === nothing ? xtimerange(first(axislimits()); tw=istw(concordiaType)) : timerange
        concordiacurve!(ax, t₀, t₁; concordiaType, kwargs...)
    end

    # Tera-Wasserburg versions, for consistency with the Plots.jl extension
    Isoplot.twconcordialine!(args...; kwargs...) = concordialine!(args...; kwargs..., concordiaType=:TeraWasserburg)
    Isoplot.twconcordiacurve!(args...; kwargs...) = concordiacurve!(args...; kwargs..., concordiaType=:TeraWasserburg)

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
