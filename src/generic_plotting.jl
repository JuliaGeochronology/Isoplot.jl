## --- Shared helpers for concordia plotting, independent of the plotting library

# Convert Wetherill (²⁰⁷Pb/²³⁵U, ²⁰⁶Pb/²³⁸U) ratios to Tera-Wasserburg (²³⁸U/²⁰⁶Pb, ²⁰⁷Pb/²⁰⁶Pb) ratios.
# Dotted so as to work equally well on scalars and arrays
twratios(r75, r68) = (1 ./ r68, r75 ./ r68 ./ value(U238_U235))

# Coordinates (in Wetherill space, or in Tera-Wasserburg space if `tw=true`) of the concordia at time(s) `t`.
# Works on numbers or arrays of ages, and on decay constants with or without uncertainties
function concordiacoords(t, λ235=λ235U, λ238=λ238U; tw::Bool=false)
    r75, r68 = ratio.(t, λ235), ratio.(t, λ238)
    return tw ? twratios(r75, r68) : (r75, r68)
end

# Intercept and slope of the straight line between times `t₀` and `t₁` on the concordia, 
# along with the x coordinates of its endpoints. (Lines in Wetherill space are also lines in Tera-Wasserburg space)
function concordialinecoeffs(t₀, t₁, λ235=λ235U, λ238=λ238U; tw::Bool=false)
    x₀, y₀ = concordiacoords(t₀, λ235, λ238; tw)
    x₁, y₁ = concordiacoords(t₁, λ235, λ238; tw)
    slope = @. (y₁ - y₀) / (x₁ - x₀)
    intercept = @. y₀ - x₀ * slope
    return intercept, slope, x₀, x₁
end

# Round-number age ticks between `tmin` and `tmax`, scaled to the order of magnitude of the time range `tlim`
function ageticks(tlim, tmin=first(tlim), tmax=last(tlim); step=0.5)
    scale = 10.0^floor(log10(tlim[2] - tlim[1])) # Characteristic timescale (order of magnitude)
    trange = round.(tlim ./ scale) # Minimum and maximum time, to a round number
    return filter(τ -> tmin <= τ <= tmax, collect((trange[1]:step:trange[2]) .* scale))
end

# Time range corresponding to the x limits `xl` of a plot, never exceeding `tmax` (by default the age of the Earth)
function xtimerange(xl; tw::Bool=false, tmax=t🜨)
    if tw # Older ages are at lower x, with x → ∞ as t → 0
        xmin = 1/ratio(tmax, value(λ238U))
        xhi = max(last(xl), xmin)
        xlo = max(first(xl), xmin)
        return age(1/xhi, value(λ238U)), age(1/xlo, value(λ238U))
    else
        t₀, t₁ = age.(max.(xl, 0.0), value(λ235U_jaffey))
        return min(t₀, tmax), min(t₁, tmax)
    end
end

# Everything needed to draw the concordia curve between times `tlim` (padded by 10% on either end):
# the curve itself, its 2σ uncertainty envelope (as a closed polygon), the same curve using the
# Schoene et al. decay constants, and age ticks, all as `(x, y)` tuples in Wetherill space or Tera-Wasserburg space.
function concordiacurvedata(tlim; tw::Bool=false, npoints::Integer=1000)
    λ235, λ238 = value(λ235U_jaffey), value(λ238U)
    # Uncertainty of 235 decay constant relative to the 238 decay constant
    σλ = λ235 * sqrt((stdev(λ238U)/λ238)^2 + (stdev(λ235U_jaffey)/λ235)^2) # 1/Myr
    dt = tlim[2] - tlim[1]
    tmin = max(tlim[1] - 0.1dt, tw ? 1e-3dt : 0.0) # Tera-Wasserburg space diverges at t=0
    tmax = tlim[2] + 0.1dt
    t = range(tmin, tmax, length=npoints)

    x, y = concordiacoords(t, λ235, λ238; tw)
    xlo, ylo = concordiacoords(t, λ235 - 2σλ, λ238; tw)
    xhi, yhi = concordiacoords(t, λ235 + 2σλ, λ238; tw)
    envelope = tw ? ([x; reverse(x)], [ylo; reverse(yhi)]) : ([xlo; reverse(xhi)], [y; reverse(y)])
    tticks = ageticks(tlim, tmin, tmax)
    return (
        curve = (x, y),
        envelope = envelope,
        schoene = concordiacoords(t, value(λ235U), λ238; tw),
        ticks = concordiacoords(tticks, λ235, λ238; tw),
        tticks = tticks,
    )
end

## --- Convert U-Pb analyses to Tera-Wasserburg space

"""
```julia
twanalysis(d::UPbAnalysis)
```
Convert a U-Pb analysis from Wetherill space (²⁰⁷Pb/²³⁵U, ²⁰⁶Pb/²³⁸U) to
Tera-Wasserburg space (²³⁸U/²⁰⁶Pb, ²⁰⁷Pb/²⁰⁶Pb), returning a two-dimensional `Analysis`.
Uncertainties and correlations are propagated to first order, and any uncertainty in 
the present-day ²³⁸U/²³⁵U ratio is ignored. The result can be plotted like any other analysis,
e.g. `plot(twanalysis.(analyses))`, and works with `datalimits` etc.

### Examples
```julia
julia> d = UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183);

julia> twanalysis(d)
Analysis2D{Float64}([2.447141738449491, 0.40236196448659206], [0.0010180454569704716, 0.00019504546431426399], [1.0364165524582163e-6 -9.273485527198267e-8; -9.273485527198267e-8 3.804273314956683e-8])
```
"""
function twanalysis(d::UPbAnalysis{T}) where {T}
    x, y = twratios(mean(d)...)
    U = T(value(U238_U235))
    # Jacobian of the transformation, column-major: x = 1/r68 and y = r75*x/U 
    J = SMatrix{2,2,T}(0, x/U, -x^2, -x*y)
    Σ = J * cov(d) * J'
    return Analysis(SVector{2,T}(x, y), (Σ + Σ')/2) # Enforce (numerically) exact symmetry
end
twanalysis(d::UThPbAnalysis) = twanalysis(UPbAnalysis(d))

## --- Weighted mean summary for rank-order plots

# Weighted mean of `data`, with the half-width of its 95% (2σ) interval, MSWD, and a legend label
function wmeanband(data::AbstractVector{<:Measurement})
    μ, mswd = wmean(data, corrected=true)
    σ = 2stdev(μ)
    digits = 1 - floor(Int, log10(σ)) # Two significant digits of uncertainty
    label = "Mean = $(round(value(μ), digits=digits)) ± $(round(σ, digits=digits)) (2σ), MSWD = $(round(mswd, sigdigits=2)), n = $(length(data))"
    return (; μ=value(μ), σ, mswd, label)
end

## --- Uncertainty ellipsoids for 3D analyses

# Surface of the covariance ellipsoid of a 3D analysis, as three `npoints`×`npoints` matrices of x, y, and z.
# Default `sigmalevel` is the trivariate p=0.05 level: sqrt(invlogccdf(Chisq(3), log(0.05)))
function ellipsoid(d::Analysis3D{T}; sigmalevel::Number=2.7954, npoints::Integer=24) where {T}
    F = eigen(Symmetric(d.Σ))
    A = F.vectors * Diagonal(T(sigmalevel) .* sqrt.(abs.(F.values))) # Maps the unit sphere onto the ellipsoid
    θ = range(0, π, length=npoints)
    φ = range(0, 2π, length=npoints)
    x, y, z = ntuple(k -> [d.μ[k] + A[k,1]*sin(t)*cos(p) + A[k,2]*sin(t)*sin(p) + A[k,3]*cos(t) for t in θ, p in φ], 3)
    return x, y, z
end

## --- Placeholder methods for plotting functions; see ext/ for (library-dependent) implementations

function concordiacurve end

function concordiacurve! end

function concordialine end

function concordialine! end

function twconcordiacurve end

function twconcordiacurve! end

function twconcordialine end

function twconcordialine! end

function rankorder end

function rankorder! end