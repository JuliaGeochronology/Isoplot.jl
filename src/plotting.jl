Plots.Shape(e::Ellipse{T}) where {T} = Shape{T,T}(e.x, e.y)

function concordiacurve!(hdl::Plots.Plot=Plots.current())
    # Uncertainty of 235 decay constant relative to the 238 decay constant
    σₜ = λ235U_jaffey.val .* sqrt((λ238U.err/λ238U.val).^2 + (λ235U_jaffey.err/λ235U_jaffey.val).^2) # 1/Years

    # Plot the concordia curve
    xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
    tlim = log.((xl .|> x-> max(x, 0.0)) .+ 1)./λ235U_jaffey.val # Calculate time range of current window
    t = range(0.9*tlim[1],1.1*tlim[2],length=1000) # Time vector, including padding
    r75t = exp.(λ235U_jaffey.val.*t) .- 1 # X axis values
    r68t = exp.(λ238U.val.*t) .- 1 # Y axis values
    x = [exp.((λ235U_jaffey.val-σₜ*2).*t) .- 1; reverse(exp.((λ235U_jaffey.val+σₜ*2).*t) .- 1)]
    y = [r68t; reverse(r68t)]
    Plots.plot!(hdl, Shape(x,y), color=:black, alpha=0.15, label="") # Two-sigma concordia uncertainty
    Plots.plot!(hdl, r75t, r68t, color=:black, label="") # Concordia line

    r75t_Schoene = exp.(λ235U.val.*t) .- 1 # X axis values
    Plots.plot!(hdl, r75t_Schoene,r68t,color=:black,style=:dash, label="") # Concordia line

    Plots.xlims!(hdl, xl) # Ensure that figure size hasn't changed
    Plots.ylims!(hdl, yl)

    # Calculate desired range of age markers
    scale=floor(log10(tlim[2]-tlim[1])) # Characteristic timescale (order of magnitude)
    trange = round.(tlim./10.0^scale) # Minimum and maximum time to a round number
    majorstep = 0.5
    tticks = (trange[1]:majorstep:trange[2]).*10.0^scale # Time ticks, to a round number
    r75tticks = exp.(λ235U_jaffey.val.*tticks) .- 1 # X axis values
    r68tticks = exp.(λ238U.val.*tticks) .- 1 # Y axis values

    # Plot age markers with text labels
    Plots.plot!(hdl, r75tticks,r68tticks, color=:black, seriestype=:scatter, ms=2, label="")
    xoffset = (xl[2]-xl[1])/200
    yoffset = (yl[2]-yl[1])/100
    t = (xl[1]+8*xoffset) .< r75tticks .< xl[2]
    ticklabels = Plots.text.(string.(round.(Int, tticks[t])), 10, :right)
    Plots.annotate!(hdl, r75tticks[t].-xoffset,r68tticks[t].+yoffset,ticklabels)

    return hdl
end

concordialine(t₀, t₁; framestyle=:box, kwargs...) = concordialine!(plot(xlims=ratio.((t₀, t₁), λ235U.val)), t₀, t₁; framestyle, kwargs...)
function concordialine!(hdl::Plots.Plot, t₀, t₁; kwargs...)
    r75₀ = ratio(t₀, λ235U.val)
    r68₀ = ratio(t₀, λ238U.val)
    r75₁ = ratio(t₁, λ235U.val)
    r68₁ = ratio(t₁, λ238U.val)
    slope = (r68₁-r68₀)/(r75₁-r75₀)
    intercept = r68₀ - r75₀*slope
    xl = Plots.xlims(hdl)
    x = collect(range(xl..., length=50))
    y = intercept .+ slope .* x
    plot!(hdl, x, val.(y); ribbon=err.(y), kwargs...)
    Plots.xlims!(hdl, xl)
end

Plots.plot(e::Ellipse, args...; kwargs...) = plot!(plot(), Shape(e), args...; kwargs...)
Plots.plot!(hdl::Plots.Plot, e::Ellipse, args...; kwargs...) = plot!(hdl, Shape(e), args...; kwargs...)
Plots.plot(e::Vector{<:Ellipse}, args...; kwargs...) = plot!(plot(), Shape.(e), args...; kwargs...)
Plots.plot!(hdl::Plots.Plot, e::Vector{<:Ellipse}, args...; kwargs...) = plot!(hdl, Shape.(e), args...; kwargs...)

Plots.plot(a::Analysis, args...; kwargs...) = plot!(plot(), ellipse(a), args...; kwargs...)
Plots.plot!(hdl::Plots.Plot, a::Analysis, args...; kwargs...) = plot!(hdl, ellipse(a), args...; kwargs...)
Plots.plot(a::Vector{<:Analysis}, args...; kwargs...) = plot!(plot(), ellipse.(a), args...; kwargs...)
Plots.plot!(hdl::Plots.Plot, a::Vector{<:Analysis}, args...; kwargs...) = plot!(hdl, ellipse.(a), args...; kwargs...)
