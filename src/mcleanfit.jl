## --- Maximum-likelihood straight-line regression in N dimensions, after McLean (2014)

# McLean, N.M. (2014) "Straight line regression through data with correlated uncertainties in
# two or more dimensions, with an application to kinetic isotope fractionation",
# Geochimica et Cosmochimica Acta 124, 237-249. doi: 10.1016/j.gca.2013.08.035

# A line is defined by a point `a` on it and a direction `v`. Given points `pᵢ` with covariance
# matrices `Σᵢ` (and `Wᵢ = Σᵢ⁻¹`), the log-likelihood of the line, profiled over the unknown position
# of each point along it, is (McLean eq. 17, dropping constants)
#     L(a, v) = -½ Σᵢ [ rᵢ'Wᵢrᵢ - (v'Wᵢrᵢ)² / (v'Wᵢv) ],   rᵢ = pᵢ - a
# The same line has many (a, v), so we fix the first component of `v` to 1 (the first variable is
# the "x" axis) and the first component of `a` to `x₀`, leaving `2(N-1)` free parameters, `θ`.
# We maximize `L` with Levenberg-Marquardt, using analytic first and second derivatives (McLean
# eqs. 18-19 and A.1-A.4). The covariance of `θ` is the inverse of the Fisher information matrix,
# i.e. of the negative Hessian at the optimum (McLean eqs. 20-21).

# In two dimensions this is the same objective as the York (1968) fit.

# The implementation here is adapted from the Matlab reference implementation by Noah McLean
# but with an interface more similar to the York fit implementation in this package,
# and internals that use SMatrix/SVector (as in `Analysis`) to reduce allocations

"""
```julia
struct McLeanFit{T,N,S}
```
Result of a maximum-likelihood straight-line fit through data in `N` dimensions, as returned by
`mcleanfit`. The fitted line is `a + v*t`, where the first component of `v` is fixed to 1 (so that
the first variable plays the role of `x`), and has fields
```
a    :: SVector{N,T}   # A point on the line, with first component equal to the weighted mean of the first variable
v    :: SVector{N,T}   # Direction of the line, with first component equal to one
Σ    :: SMatrix        # Covariance of the 2(N-1) free components of `a` and `v` (excluding the first of each)
mswd :: T              # Reduced chi-squared, with (N-1)(n-2) degrees of freedom
```
For convenience, the properties `slope` and `intercept` (the free components of `v`, and the
point where the first variable is zero) give 1σ `Measurement`s, scalar if `N == 2` and
`SVector`s otherwise. In two dimensions the properties `xm` and `ym` are also available, giving a
point on the line at which the uncertainties in `ym` and `slope` are uncorrelated, just as for
`YorkFit`.
"""
struct McLeanFit{T<:Number, N, S<:SMatrix} <: Data{T}
    a::SVector{N,T}
    v::SVector{N,T}
    Σ::S
    mswd::T
end
# Other constructors
McLeanFit(::Val{N}, x::F) where {N,F} = McLeanFit(
    fill(x, SVector{N,F}),
    fill(x, SVector{N,F}),
    fill(x, SMatrix{2N-2,2N-2,F,(2N-2)^2}),
    x,
)

# Derived properties
Base.propertynames(::McLeanFit{<:Any,2}) = (:a, :v, :Σ, :mswd, :intercept, :slope, :xm, :ym)
Base.propertynames(::McLeanFit) = (:a, :v, :Σ, :mswd, :intercept, :slope)
function Base.getproperty(f::McLeanFit{T,N}, s::Symbol) where {T,N}
    if s === :slope || s === :intercept
        # Marginal uncertainties of the free components, which are the first `N-1` (a) and last `N-1` (v) of Σ
        a, v, Σ = getfield(f, :a), getfield(f, :v), getfield(f, :Σ)
        ia = SVector{N-1}(ntuple(identity, Val(N-1)))
        iv = ia .+ (N-1)
        vᵥ = v[ia .+ 1]
        if s === :slope
            m = vᵥ .± sqrt.(diag(Σ[iv,iv]))
        else # intercept, i.e. where first component is zero: aⱼ - a₁vⱼ
            aᵥ = a[ia .+ 1] .- a[1] .* vᵥ
            σ² = diag(Σ[ia,ia]) .+ a[1]^2 .* diag(Σ[iv,iv]) .- 2a[1] .* diag(Σ[ia,iv])
            m = aᵥ .± sqrt.(σ²)
        end
        return N == 2 ? only(m) : m
    elseif s === :xm || s === :ym
        N == 2 || throw(ArgumentError("`$s` is only defined for two-dimensional fits"))
        a, b, Σ = getfield(f, :a), getfield(f, :v)[2], getfield(f, :Σ)
        # Shift reference point along the line so that uncertainties in its height and the slope are uncorrelated
        δ = -Σ[1,2] / Σ[2,2]
        s === :xm && return a[1] + δ
        return (a[2] + b*δ) ± sqrt(Σ[1,1] - Σ[1,2]^2 / Σ[2,2])
    else
        return getfield(f, s)
    end
end

# Calculate a line given data (see also `line` and `invline` for `YorkFit`s)
line(f::McLeanFit{<:Any,2}, x::Number) = f.ym + (x - f.xm) * f.slope
invline(f::McLeanFit{<:Any,2}, y::Number) = (y - f.ym) / f.slope + f.xm

"""
```julia
mcleanfit(x::Collection, σx::Collection, y::Collection, σy::Collection, [ρ]; maxiter=128, tol=1e-12, corrected=false)
mcleanfit(x::Collection{<:Measurement}, y::Collection{<:Measurement}, [ρ]; kwargs...)
mcleanfit(d::Collection{<:AbstractAnalysis}; kwargs...)
mcleanfit(μ::Collection{<:AbstractVector}, Σ::Collection{<:AbstractMatrix}; kwargs...)
```
Fit a straight line through data in two or more dimensions, with arbitrary correlated
uncertainties, following the maximum-likelihood Levenberg-Marquardt approach of McLean (2014),
and return a `McLeanFit`. The first variable is treated as the independent (`x`) variable.

Data may be given as in `yorkfit` for two dimensions: `x`, `y`, one-sigma uncertainties `σx` and
`σy`, and optionally correlation coefficients `ρ`. Alternatively, in any number of dimensions,
data may be given as a collection of `Analysis` objects (or of `UPbAnalysis` etc.), or as
collections of means `μ` and covariance matrices `Σ`. Points with any `NaN`s are ignored.

In two dimensions the fitted line is equivalent to that of `yorkfit`.

Uncertainties are 1σ. If `corrected=true`, they are expanded by `sqrt(mswd)` whenever `mswd > 1`.
The solver stops when the expected improvement in chi-squared from a further Newton step falls
below `tol`, or after `maxiter` iterations, whichever happens first.

For further reference, see:
McLean, N.M. (2014) "Straight line regression through data with correlated uncertainties in
two or more dimensions, with an application to kinetic isotope fractionation"
Geochimica et Cosmochimica Acta 124, 237-249. doi: 10.1016/j.gca.2013.08.035

### Examples
```julia
julia> x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588];

julia> y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955];

julia> mcleanfit(x, fill(0.25, 10), y, fill(0.25, 10))
McLeanFit{Float64}:
 Maximum-likelihood linear fit of the form y = a + bx with
  intercept: -0.24 ± 0.25 (1σ)
  slope    : 1.043 ± 0.04 (1σ)
  MSWD     : 1.4272439423337353
```
"""
function mcleanfit(x::Collection{T}, σx::Collection, y::Collection, σy::Collection, ρ=zero(T); kwargs...) where {T}
    @assert eachindex(x) == eachindex(σx) == eachindex(y) == eachindex(σy) "Indices of x, σx, y, and σy must match"
    @assert !isempty(x) "Mclean Fit requires at least one data point"
    @assert ρ isa Number || eachindex(ρ) == eachindex(x) "Indices of correlation coefficients must match data"
    @assert all(x->(-1 <= x <= 1), ρ) "Correlation coefficients must be between -1 and 1"
    F = float(T)
    p = [SVector{2,F}(x[i], y[i]) for i in eachindex(x)]
    Σ = map(eachindex(x)) do i
        ρᵢ = _maybe_index(ρ, i)
        SMatrix{2,2,F,4}(σx[i]^2, ρᵢ*σx[i]*σy[i], ρᵢ*σx[i]*σy[i], σy[i]^2)
    end
    return mcleanfit(p, Σ; kwargs...)
end
mcleanfit(x::Collection{<:Measurement}, y::Collection{<:Measurement}, ρ=ntuple(i->paircor(x[i], y[i]), length(x)); kwargs...) = mcleanfit(value.(x), stdev.(x), value.(y), stdev.(y), ρ; kwargs...)
function mcleanfit(d::Collection{<:AbstractAnalysis}; kwargs...)
    @assert !isempty(d)
    return mcleanfit([mean(dᵢ) for dᵢ in d], [cov(dᵢ) for dᵢ in d]; kwargs...)
end
mcleanfit(p::Collection{<:SVector}, Σ::Collection{<:SMatrix}; kwargs...) = mcleanfit(collect(p), collect(Σ); kwargs...)
function mcleanfit(μ::Collection{<:AbstractVector}, Σ::Collection{<:AbstractMatrix}; kwargs...)
    N, F = length(first(μ)), float(promote_type(eltype(first(μ)), eltype(first(Σ))))
    return mcleanfit([SVector{N,F}(m) for m in μ], [SMatrix{N,N,F}(C) for C in Σ]; kwargs...)
end

# Core method, given vectors of points and covariance matrices
function mcleanfit(p::Vector{SVector{N,F}}, Σ::Vector{SMatrix{N,N,F,N²}}; maxiter::Integer=128, tol::Number=1e-12, corrected::Bool=false) where {N,F<:AbstractFloat,N²}
    @assert N >= 2 "Data must have at least two dimensions"
    @assert eachindex(p) == eachindex(Σ)

    # Ignore any points with NaNs
    isbad(i) = any(isnan, p[i]) || any(isnan, Σ[i])
    if any(isbad, eachindex(p))
        ok = [!isbad(i) for i in eachindex(p)]
        return mcleanfit(p[ok], Σ[ok]; maxiter, tol, corrected)
    end
    n = length(p)
    (n < 2 || !all(isposdef, Σ)) && return McLeanFit(Val(N), F(NaN))
    W = map(inv, Σ)

    # Fix first component of `a` to its value in the multivariate weighted mean (McLean, appendix B)
    p̂ = sum(W) \ sum(W[i]*p[i] for i in eachindex(p))
    x₀ = p̂[1]

    # Initial estimate of the line from OLS regression of each variable against the first one
    p̄ = sum(p) / n
    Sxx = sum(q->(q[1]-p̄[1])^2, p)
    v₀ = sum(q->(q[1]-p̄[1])*(q-p̄), p) / Sxx # First component is exactly one
    (isfinite(Sxx) && Sxx > 0) || return McLeanFit(Val(N), F(NaN))
    a₀ = p̄ + v₀*(x₀ - p̄[1])

    # Maximum likelihood estimate of the line parameters
    S, z₀ = _mclean_parameterization(x₀, Val(N))
    θ, L, H = mclean_solve(S'*vcat(a₀, v₀), S, z₀, p, W, F(tol), maxiter)
    isnan(L) && return McLeanFit(Val(N), F(NaN))
    isposdef(Symmetric(-H)) || return McLeanFit(Val(N), F(NaN))

    # Covariance of the fitted parameters, expanded for excess dispersion if requested
    mswd = n > 2 ? max(-2L / ((N-1)*(n-2)), zero(F)) : F(NaN)
    C = inv(-H)
    C = (C + C') / 2 * (corrected ? max(mswd, one(F)) : one(F))
    a, v = _mclean_halves(z₀ + S*θ, Val(N))
    return McLeanFit(a, v, C, mswd)
end

# Free parameters are all components of `a` and `v` but the first of each. Returns the
# 2N × 2(N-1) matrix `S` and the length 2N vector `z₀` such that `z = (a, v) = z₀ + S*θ`
function _mclean_parameterization(x₀::F, ::Val{N}) where {F,N}
    S = SMatrix{2N,2N-2,F}(ntuple(Val(2N*(2N-2))) do k
        i, j = (k-1) % 2N + 1, (k-1) ÷ 2N + 1
        F(i == j + 1 + (j >= N)) # a₂...a_N map to rows 2:N, v₂...v_N to rows N+2:2N
    end)
    z₀ = SVector{2N,F}(ntuple(k -> k == 1 ? x₀ : (k == N+1 ? one(F) : zero(F)), Val(2N)))
    return S, z₀
end
_mclean_halves(z::SVector{M,F}, ::Val{N}) where {M,F,N} = SVector{N,F}(ntuple(i->z[i], Val(N))), SVector{N,F}(ntuple(i->z[N+i], Val(N)))

# Levenberg-Marquardt iterations, starting with little damping. Returns the optimum `θ`, the
# log-likelihood `L` there, and the Hessian `H` of `L` with respect to `θ`.
function mclean_solve(θ::SVector{M,F}, S::SMatrix{M2,M,F}, z₀, p, W, tol, maxiter) where {M,M2,F}
    N = Val(M2÷2)
    Id = SMatrix{M,M,F}(I)
    # We're done once the expected gain in log-likelihood from a full Newton step (the Newton
    # decrement, in units of χ²) is negligible. Since the χ² of each point is a difference of
    # terms that can be much larger than χ² itself, we can't expect to resolve it below roundoff
    # error in those terms
    a₀, _ = _mclean_halves(z₀ + S*θ, N)
    scale = sum(i -> (ρᵢ = p[i] - a₀; ρᵢ ⋅ (W[i]*ρᵢ)), eachindex(p))
    tolₑ = max(tol, 100*eps(F)*scale)

    function evaluate(θ, derivatives=true)
        a, v = _mclean_halves(z₀ + S*θ, N)
        derivatives || return _mclean_loglikelihood(a, v, p, W)
        Gz, Hz = _mclean_derivatives(a, v, p, W)
        return S' * Gz, S' * Hz * S # Restrict to the free parameters
    end
    λ = F(1e-3)
    L = evaluate(θ, false)
    G, H = evaluate(θ)
    converged = false
    for _ in 1:maxiter
        # Only a maximum, with negative definite Hessian, counts as converged (not a saddle point)
        if isposdef(Symmetric(-H)) && -G ⋅ (H\G) / 2 <= tolₑ
            converged = true
            break
        end
        # Marquardt damping, increased until the damped Hessian is negative definite so that steps go uphill
        Hₗ = H - λ*(abs.(H) .* Id)
        if !isposdef(Symmetric(-Hₗ))
            λ *= 10
            isfinite(λ) || break
            continue
        end
        h = -Hₗ \ G
        Lₙ = evaluate(θ + h, false)
        if Lₙ > L
            # Accept step and relax damping
            θ += h
            L = Lₙ
            λ /= 10
            G, H = evaluate(θ)
        else
            # Reject step and increase damping (also catches NaNs)
            λ *= 10
            isfinite(λ) || break
        end
    end
    converged || @warn "McLean fit failed to converge after $maxiter iterations"

    # Polish with one last Newton step if converged to tolerance (cheap, and quadratically convergent from here)
    converged && (θ -= H \ G)
    return θ, evaluate(θ, false), last(evaluate(θ))
end

# Profile log-likelihood (minus constants) of the line `a + vt`, McLean eq. 17
function _mclean_loglikelihood(a::SVector{N,F}, v::SVector{N,F}, p, W) where {N,F}
    χ² = zero(F)
    @inbounds for i in eachindex(p)
        ρᵢ = p[i] - a
        t, u = W[i]*ρᵢ, W[i]*v
        χ² += ρᵢ ⋅ t - (v ⋅ t)^2 / (v ⋅ u)
    end
    return -χ²/2
end

# Gradient and Hessian of the above with respect to all of `(a, v)`, McLean eqs. 18, 19, and A.1-A.4
function _mclean_derivatives(a::SVector{N,F}, v::SVector{N,F}, p, W) where {N,F}
    ∂a, ∂v = zeros(SVector{N,F}), zeros(SVector{N,F})
    ∂a∂a, ∂a∂v, ∂v∂v = zeros(SMatrix{N,N,F}), zeros(SMatrix{N,N,F}), zeros(SMatrix{N,N,F})
    @inbounds for i in eachindex(p,W)
        Wᵢ = W[i]
        rᵢ = p[i] - a
        t, u = Wᵢ*rᵢ, Wᵢ*v # Terms that appear throughout
        q, s = v ⋅ u, v ⋅ t
        P = Wᵢ - u*u'/q
        ∂a += P*rᵢ
        ∂v += (s/q)*t - (s/q)^2*u
        ∂a∂a -= P
        ∂a∂v -= (q*(u*t' + s*Wᵢ) - 2s*u*u') / q^2
        ∂v∂v += (q^2*(2s*t*u' + q*t*t' - s^2*Wᵢ - 2s*u*t') - 4q*(q*s*t - s^2*u)*u') / q^4
    end
    return vcat(∂a, ∂v), vcat(hcat(∂a∂a, ∂a∂v), hcat(∂a∂v', ∂v∂v))
end
