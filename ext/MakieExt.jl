module MakieExt

    using Isoplot
    using Makie
    using Measurements
    import Isoplot: concordialine, concordialine!, concordiacurve, concordiacurve!
    import Isoplot: concordiacoords, concordialinecoeffs, concordiacurvedata, xtimerange, wmeanband, ellipsoid

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

    # 3D analyses are drawn as translucent ellipsoids; use `Ellipse(a, i, j)` for a 2D projection onto variables i and j
    @recipe(Ellipsoids, a) do scene
        Attributes(
            color = Makie.wong_colors()[1],
            alpha = 0.3,
            sigmalevel = 2.7954,
            npoints = 24,
            label = nothing,
        )
    end
    Makie.convert_arguments(::Type{<:Ellipsoids}, a::Union{Analysis3D, AbstractVector{<:Analysis3D}}) = (a,)
    Makie.plottype(::Analysis3D) = Ellipsoids
    Makie.plottype(::AbstractVector{<:Analysis3D}) = Ellipsoids
    Makie.legendelements(p::Ellipsoids, legend) = [PolyElement(color=(p.color[], p.alpha[]))]

    function Makie.plot!(p::Ellipsoids)
        as = p.a[] isa Analysis3D ? [p.a[]] : p.a[]
        for a in as
            x, y, z = ellipsoid(a; sigmalevel=p.sigmalevel[], npoints=p.npoints[])
            surface!(p, x, y, z; color=fill(p.color[], size(x)), alpha=p.alpha[], transparency=true, shading=NoShading)
        end
        return p
    end

    function as_points(v1, v2)
        length(v1) == length(v2) || throw(DimensionMismatch("Arguments must be the same length"))
        return Point2f[Point2f(v1[i], v2[i]) for i in 1:length(v1)] # (Not `zip`, which would make a statically-sized result for SVectors)
    end

    ## --- Helpers

    # Current limits of the active axis, as a tuple of `(min, max)` for each dimension
    function axislimits()
        ax = Makie.current_axis()
        ax isa Union{Makie.Axis, Makie.Axis3} || throw(ArgumentError("No current Axis to take limits from; pass limits explicitly"))
        Makie.reset_limits!(ax)
        r = ax.finallimits[]
        return ntuple(i -> (r.origin[i], r.origin[i] + r.widths[i]), length(r.origin))
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

    ## --- Line fits (York and McLean)

    # The axis decides the dimensionality, as in Makie generally. On an `Axis`, `dims` is an integer: plot that
    # variable against the first, with 2σ uncertainty. On an `Axis3`, `dims` is a tuple of two integers: plot a
    # line through those variables against the first (no uncertainty).
    @recipe(LineFit, lf) do scene
        Attributes(
            color = Makie.wong_colors()[2],
            linewidth = 2,
            bandalpha = 0.3,
            npoints = 25,
            dims = nothing, # Defaults to 2 on an `Axis`, or (2,3) on an `Axis3`
            xlims = nothing, # Defaults to the limits of the current axis
            label = nothing,
        )
    end
    Makie.convert_arguments(::Type{<:LineFit}, lf::Union{YorkFit, McLeanFit}) = (lf,)
    Makie.plottype(::Union{YorkFit, McLeanFit}) = LineFit

    function Makie.plot!(p::LineFit)
        lf = p.lf[]
        dims = p.dims[] === nothing ? (Makie.current_axis() isa Makie.Axis3 ? (2,3) : 2) : p.dims[]
        xl = p.xlims[] === nothing ? first(axislimits()) : p.xlims[]
        x = collect(range(xl..., length=p.npoints[]))
        if dims isa Integer
            y = Isoplot.line.(lf, x, dims)
            ŷ, σ = value.(y), 2stdev.(y)
            band!(p, x, ŷ .- σ, ŷ .+ σ; color=(p.color[], p.bandalpha[]))
            lines!(p, x, ŷ; color=p.color, linewidth=p.linewidth, label=p.label)
        else
            y, z = (value.(Isoplot.line.(lf, x, d)) for d in dims)
            lines!(p, x, y, z; color=p.color, linewidth=p.linewidth, label=p.label)
        end
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
    # With `wmean=true`, also overlay the weighted mean and its 2σ interval
    function Isoplot.rankorder!(ax::Makie.AbstractAxis, args...; wmean=false, meancolor=:black, meanlabel=:auto, scale=1, offset=0, kwargs...)
        p = rankorderplot!(ax, args...; scale, offset, kwargs...)
        if wmean
            data = p.data[]
            data isa AbstractVector{<:Measurement} || throw(ArgumentError("`wmean` is only supported for `Measurement` data"))
            wmeanplot!(ax, data; scale, offset, color=meancolor, label=meanlabel === :auto ? wmeanband(data).label : meanlabel)
        end
        return p
    end

    # Weighted mean of rank-ordered data, as a line with a 2σ band
    @recipe(WMeanPlot, data) do scene
        Attributes(color = :black, scale = 1, offset = 0, label = nothing)
    end
    Makie.legendelements(p::WMeanPlot, legend) = [PolyElement(color=(p.color[], 0.15)), LineElement(color=p.color[])]

    function Makie.plot!(p::WMeanPlot)
        data = p.data[]
        m = wmeanband(data)
        n = length(data)
        xs = [p.offset[] + p.scale[]*(1 - 1/2), p.offset[] + p.scale[]*(n + 1/2)]
        band!(p, xs, fill(m.μ - m.σ, 2), fill(m.μ + m.σ, 2); color=(p.color[], 0.15))
        lines!(p, xs, fill(m.μ, 2); color=p.color, label=p.label)
        return p
    end

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
