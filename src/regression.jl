using SpecialFunctions: erfc

"""
```julia
chauvenet(μ::AbstractVector, σ::AbstractVector)
```
Apply Chauvenet's criterion to a set of data to exclude outliers.

The function calculates the z-scores of the data points, and then calculates the probability `p` of observing a value as extreme as the z-score under the assumption of normal distribution.
It then applies Chauvenet's criterion, marking any data point as an outlier if `2 * N * p < 1.0`, where `N` is the total number of data points.

### Examples
```julia

```
"""
function chauvenet(μ::AbstractVector, σ::AbstractVector; threshold=1.0)
    @assert eachindex(μ) == eachindex(σ)
    criterion = chauvenet_criterion(μ, σ)
    selected = criterion .>= threshold
    return μ[selected], σ[selected]
end
function chauvenet(x::AbstractVector; threshold=1.0)
    criterion = chauvenet_criterion(value.(x), stdev.(x))
    return x[criterion .>= threshold]
end
function chauvenet_criterion(μ::AbstractVector, σ::AbstractVector)
    mean_val = mean(μ)
    N = length(μ)
    z_scores = abs.(μ .- mean_val) ./ σ
    p = 0.5 * erfc.(z_scores ./ sqrt(2.0))
    criterion = 2 * N * p
    return criterion
end

## --- Weighted means
"""
```julia
wμ, wσ, mswd = wmean(μ, σ; corrected=true)
wμ ± wσ, mswd = wmean(μ ± σ; corrected=true)
```
The weighted mean, with or without the "geochronologist's MSWD correction" to uncertainty.
You may specify your means and standard deviations either as separate vectors `μ` and `σ`,
or as a single vector `x` of `Measurement`s equivalent to `x = μ .± σ`

In all cases, `σ` is assumed to reported as _actual_ sigma (i.e., 1-sigma).

If `corrected=true`, the resulting uncertainty of the weighted mean is expanded by a factor
of `sqrt(mswd)` to attempt to account for excess dispersion when the MSWD is greater than `1`

### Examples
```julia
julia> x = randn(10)
10-element Vector{Float64}:
  0.4612989881720301
 -0.7255529837975242
 -0.18473979056481055
 -0.4176427262202118
 -0.21975911391551833
 -1.6250003193791873
 -1.6185557291787287
  0.25315988825847513
 -0.4979804844182867
  1.3565281078086726

julia> y = ones(10);

julia> wmean(x, y)
(-0.321824416323509, 0.31622776601683794, 0.8192171477885678)

julia> wmean(x .± y)
(-0.32 ± 0.32, 0.8192171477885679)

julia> wmean(x .± y./10)
(-0.32 ± 0.29, 81.92171477885678)

julia> wmean(x .± y./10, corrected=false)
(-0.322 ± 0.032, 81.92171477885678)
```
"""
# Weighted mean given a vector of means and vector of standard deviations, assuming no covariance
function wmean(μ::Collection1D{T}, σ::Collection1D{T}; corrected::Bool=true) where {T}
    sum_of_values = sum_of_weights = χ² = zero(float(T))
    @inbounds for i in eachindex(μ,σ)
        σ² = σ[i]^2
        sum_of_values += μ[i] / σ²
        sum_of_weights += one(T) / σ²
    end
    wμ = sum_of_values / sum_of_weights
    wσ² = 1 / sum_of_weights

    @inbounds for i in eachindex(μ,σ)
        χ² += (μ[i] - wμ)^2 / σ[i]^2
    end
    mswd = χ² / (length(μ)-1)

    # Optional: expand standard error by sqrt of mswd, if mswd > 1
    corrected && (wσ² *= nanmax(mswd,1))
    
    return wμ, sqrt(wσ²), mswd
end
# Full covariance matrix method
function wmean(x::AbstractVector{T}, C::AbstractMatrix{T}; corrected::Bool=true) where T
    # Weighted mean and variance, full matrix method
    J = ones(length(x))
    σ²ₓ̄ = 1/(J'/C*J)
    x̄ = σ²ₓ̄*(J'/C*x)

    # MSWD, full matrix method
    r = x .- x̄
    χ² = r'/C*r
    mswd = χ² / (length(x)-1)

    # Optional: expand standard error by sqrt of mswd, if mswd > 1
    corrected && (σ²ₓ̄ *= nanmax(mswd,1))

    return x̄, sqrt(σ²ₓ̄), mswd
end
# Weighted mean given a vector of Distributions, assuming no covariance
function wmean(analyses::Collection{<:Distribution}; corrected::Bool=true)
    T = eltype(eltype(analyses))
    sum_of_values = sum_of_weights = χ² = zero(T)
    @inbounds for i in eachindex(analyses)
        σ² = var(analyses[i])
        sum_of_values += mean(analyses[i]) / σ²
        sum_of_weights += one(T) / σ²
    end
    wμ = sum_of_values / sum_of_weights
    wσ² = 1 / sum_of_weights

    @inbounds for i in eachindex(analyses)
        χ² += (mean(analyses[i]) - wμ)^2 / std(analyses[i])^2
    end
    mswd = χ² / (length(analyses)-1)

    # Optional: expand standard error by sqrt of mswd, if mswd > 1
    corrected && (wσ² *= nanmax(mswd,1))

    return Normal(wμ, sqrt(wσ²)), mswd
end
# Weighted mean given a vector of Measurements (possibly with covariance)
function wmean(x::AbstractVector{Measurement{T}}; corrected::Bool=true) where {T}
    wμ, wσ, mswd = wmean(value.(x), Measurements.cov(x); corrected)
    return wμ ± wσ, mswd
end

# Legacy methods, for backwards compatibility
awmean(args...) = wmean(args...; corrected=false)
gwmean(args...) = wmean(args...; corrected=true)


distwmean(x...; corrected::Bool=true) = distwmean(x; corrected)
function distwmean(x::Collection{<:Distribution}; corrected::Bool=true, N=10^6)
    @assert eachindex(x) == 1:length(x)
    sampled = ntuple(i->rand(x[i], N), length(x))
    distwmean(sampled; corrected)
end
function distwmean(x::NTuple{N, <:AbstractVector}; corrected::Bool=true) where {N}
    σₓ = nanstd.(x)
    wₓ = 1 ./ σₓ.^2
    wₜ = sum(wₓ)
    c = similar(first(x))
    @inbounds for i in eachindex(x...)
        c[i] = sum(getindex.(x, i) .* wₓ)/wₜ
        if corrected
            mswd = sum(wₓ.*abs2.(getindex.(x,i).-c[i]))/(wₜ*(N-1))
            c[i] += sqrt(mswd) * randn()
        end
    end
    return c
end

"""
```julia
mswd(μ, σ)
mswd(μ ± σ)
```
Return the Mean Square of Weighted Deviates (AKA the reduced chi-squared
statistic) of a dataset with values `x` and one-sigma uncertainties `σ`

### Examples
```julia
julia> x = randn(10)
10-element Vector{Float64}:
 -0.977227094347237
  2.605603343967434
 -0.6869683962845955
 -1.0435377148872693
 -1.0171093080088411
  0.12776158554629713
 -0.7298235147864734
 -0.3164914095249262
 -1.44052961622873
  0.5515207382660242

julia> mswd(x, ones(10))
1.3901517474017941
```
"""
function mswd(μ::Collection1D{T}, σ::Collection) where {T}
    sum_of_values = sum_of_weights = χ² = zero(float(T))

    @inbounds for i in eachindex(μ,σ)
        w = 1 / σ[i]^2
        sum_of_values += w * μ[i]
        sum_of_weights += w
    end
    wx = sum_of_values / sum_of_weights

    @inbounds for i in eachindex(μ,σ)
        χ² += (μ[i] - wx)^2 / σ[i]^2
    end

    return χ² / (length(μ)-1)
end
# Full covariance matrix method
function mswd(x::AbstractVector{T}, C::AbstractMatrix{T}) where T
    # Weighted mean and variance, full matrix method
    J = ones(length(x))
    σ²ₓ̄ = 1/(J'/C*J)
    x̄ = σ²ₓ̄*(J'/C*x)

    # MSWD, full matrix method
    r = x .- x̄
    χ² = r'/C*r
    return χ² / (length(x)-1)
end
mswd(x::AbstractVector{<:Measurement}) = mswd(value.(x), Measurements.cov(x))
mswd(x::AbstractVector) = mswd(value.(x), stdev.(x))

## ---  Simple linear regression

"""
```julia
(a,b) = lsqfit(x::AbstractVector, y::AbstractVector)
```
Returns the coefficients for a simple linear least-squares regression of
the form `y = a + bx`

### Examples
```
julia> a, b = lsqfit(1:10, 1:10)
2-element Vector{Float64}:
 -1.19542133983862e-15
  1.0

julia> isapprox(a, 0, atol = 1e-12)
true

julia> isapprox(b, 1, atol = 1e-12)
true
```
"""
lsqfit(x::Collection{<:Number}, y::Collection{<:Number}) = lsqfit(x, collect(y))
function lsqfit(x::Collection{T}, y::AbstractVector{<:Number}) where {T<:Number}
    A = Array{T}(undef, length(x), 2)
    A[:,1] .= one(T)
    A[:,2] .= x
    return A\y
end
# Identical to the one in StatGeochemBase

## -- The York (1968) two-dimensional linear regression with x and y uncertainties
    # as commonly used in isochrons

# Custom type to hold York fit resutls
struct YorkFit{T<:Number} <: Data{T}
    intercept::Measurement{T}
    slope::Measurement{T}
    xm::T
    ym::Measurement{T}
    mswd::T
end

# Calculate a line given data
line(yf::YorkFit, x::Number) =  yf.ym + (x - yf.xm) * yf.slope
invline(yf::YorkFit, y::Number) = (y - yf.ym) / yf.slope + yf.xm

"""
```julia
yorkfit(x::Collection, σx::Collection, y::Collection, σy::Collection, [r]; maxiter=128, tol=1e-12, corrected=false)
yorkfit(x::Collection{<:Measurement}, y::Collection{<:Measurement}, [r]; maxiter=128, tol=1e-12, corrected=false)
yorkfit(d::Collection{<:AbstractAnalysis}; maxiter=128, tol=1e-12, corrected=false)
```
Uses the York (1968) two-dimensional iterative least-squares fit to calculate `a`, `b`,
and uncertanties `σa`, `σb` for the equation `y = a + bx`, given `x`, `y`,
uncertaintes `σx`, and `σy`, and optionally covarances `r`.
Iteration continues until the slope converges to relative tolerance `tol`,
or for at most `maxiter` iterations, whichever comes first.

If `corrected=true`, the resulting uncertainties in slope (`σb`), intercept (`σa`),
and midpoint (`σym`) are expanded by a factor of `sqrt(mswd)` 
to attempt to account for excess dispersion when the MSWD is greater than `1`.

For further reference, see:
York, Derek (1968) "Least squares fitting of a straight line with correlated errors"
Earth and Planetary Science Letters 5, 320-324. doi: 10.1016/S0012-821X(68)80059-7

### Examples
```julia
julia> x = (1:100) .+ randn.();

julia> y = 2*(1:100) .+ randn.();

julia> yorkfit(x, ones(100), y, ones(100))
YorkFit{Float64}:
 Least-squares linear fit of the form y = a + bx with
  intercept: -0.54 ± 0.45 (1σ)
  slope    : 2.0065 ± 0.0078 (1σ)
  MSWD     : 0.9457955207778813
```
"""
function yorkfit(d::Collection{<:AbstractAnalysis}; kwargs...)
    @assert !isempty(d)
    # Using NTuples instead of Arrays can avoid allocations if length(d) can be const-prop'd
    # and should be more efficient for relatively small N, but could be less efficient for large N
    x = ntuple(i->mean(d[i])[1], length(d))
    y = ntuple(i->mean(d[i])[2], length(d))
    σx = ntuple(i->std(d[i])[1], length(d))
    σy = ntuple(i->std(d[i])[2], length(d))
    r = ntuple(i->cov(d[i])[1,2]/(sqrt(cov(d[i])[1,1] * cov(d[i])[2,2])), length(d))
    yorkfit(x, σx, y, σy, r; kwargs...)
end
yorkfit(x::Collection{<:Measurement}, y::Collection{<:Measurement}, r=ntuple(i->paircor(x[i], y[i]), length(x)); kwargs...) = yorkfit(value.(x), stdev.(x), value.(y), stdev.(y), r; kwargs...)
function yorkfit(x::Collection{T}, σx::Collection, y::Collection{T}, σy::Collection, r=zero(T); maxiter::Integer=128, tol::Number=1e-12, corrected::Bool=false) where {T}
    @assert !isempty(x) "Cannot take York fit of empty data"
    @assert eachindex(x) == eachindex(σx) == eachindex(y) == eachindex(σy) "Indices of x, σx, y, and σy must match"
    @assert r isa Number || eachindex(r) == eachindex(x) "Indices of correlation coefficients must match data"
    @assert all(x->(-1<=x<=1), r) "Correlation coefficients must be between -1 and 1"
    n = length(x)

    ## For an initial estimate of the slope, calculate the ordinary
    # least-squares fit for the equation y=a+bx
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    Σxy, Σxx = zero(x̄*ȳ), zero(x̄*x̄)
    @inbounds for i in eachindex(x)
        Σxy += (x[i] - x̄) * (y[i] - ȳ)
        Σxx += (x[i] - x̄)^2
    end
    b = Σxy / Σxx
    ∅ = zero(b)

    ## Perform the York fit, iterating the slope to convergence
    # Each pass recomputes the weights `Z` from the current slope, and then
    # the weighted means and the updated slope.
    tolₑ = max(tol, 10*eps(typeof(b)))
    converged = false
    for _ in 1:maxiter
        ΣZ, ΣZx, ΣZy = ∅, ∅, ∅
        @inbounds for i in eachindex(x)
            Z = _york_weight(b, σx[i], σy[i], _maybe_index(r, i))
            ΣZ += Z
            ΣZx += Z * x[i]
            ΣZy += Z * y[i]
        end
        x̄ = ΣZx / ΣZ
        ȳ = ΣZy / ΣZ

        sV, sU = ∅, ∅
        @inbounds for i in eachindex(x)
            ωx, ωy = inv(σx[i]^2), inv(σy[i]^2)
            ρ = _maybe_index(r, i)
            α = sqrt(ωx * ωy)
            Z = ωx*ωy / (b^2*ωy + ωx - 2*b*ρ*α)
            U, V = x[i] - x̄, y[i] - ȳ
            sV += Z^2 * V * (U/ωy + b*V/ωx - ρ*V/α)
            sU += Z^2 * U * (U/ωy + b*V/ωx - b*ρ*U/α)
        end
        bₙ = sV / sU
        converged = abs(bₙ - b) <= tolₑ * abs(bₙ)
        b = bₙ
        converged && break
    end
    converged || @warn "York fit did not converge within $maxiter iterations; consider `mcleanfit`, which may be more robust"

    ## Weighted means and Z at the final slope
    ΣZ, ΣZx, ΣZy = ∅, ∅, ∅
    @inbounds for i in eachindex(x)
        Z = _york_weight(b, σx[i], σy[i], _maybe_index(r, i))
        ΣZ += Z
        ΣZx += Z * x[i]
        ΣZy += Z * y[i]
    end
    x̄ = ΣZx / ΣZ
    ȳ = ΣZy / ΣZ
    a = ȳ - b * x̄

    ## Adjusted points (u, v) and their weighted mean (xm, ym)
    ΣZu, Σχ² = ∅, ∅
    @inbounds for i in eachindex(x)
        Z, u = _york_adjusted(b, x̄, ȳ, x[i], y[i], σx[i], σy[i], _maybe_index(r, i))
        ΣZu += Z * u
        Σχ² += (y[i] - a - b*x[i])^2 / (σy[i]^2 + b^2*σx[i]^2 - 2*b*_maybe_index(r, i)*σx[i]*σy[i])
    end
    xm = ΣZu / ΣZ
    ym = ∅
    @inbounds for i in eachindex(x)
        Z, u = _york_adjusted(b, x̄, ȳ, x[i], y[i], σx[i], σy[i], _maybe_index(r, i))
        ym += Z * (ȳ + b*(u - x̄))
    end
    ym /= ΣZ

    ## Uncertainties and MSWD
    Σd = ∅
    @inbounds for i in eachindex(x)
        Z, u = _york_adjusted(b, x̄, ȳ, x[i], y[i], σx[i], σy[i], _maybe_index(r, i))
        Σd += Z * (u - xm)^2
    end
    σb = sqrt(inv(Σd))
    σa = sqrt(inv(ΣZ) + xm^2 * σb^2)
    σym = sqrt(inv(ΣZ))

    # MSWD (reduced chi-squared) of the fit
    mswd = Σχ² / (n-2)

    # Optionally scale uncertainties by √MSWD if the data are overdispersed
    if corrected
        f = sqrt(max(mswd, one(mswd)))
        σb *= f
        σa *= f
        σym *= f
    end

    ## Results
    return YorkFit(a ± σa, b ± σb, xm, ym ± σym, mswd)
end

# Helpers for `yorkfit`
@inline _maybe_index(r::Number, i) = r
@inline _maybe_index(r, i) = @inbounds r[i]
@inline function _york_weight(b, σx, σy, r)
    ωx, ωy = inv(σx^2), inv(σy^2)
    return ωx*ωy / (b^2*ωy + ωx - 2*b*r*sqrt(ωx*ωy))
end
# Weight `Z` and adjusted x-position `u` of a point, given the slope and weighted means
@inline function _york_adjusted(b, x̄, ȳ, x, y, σx, σy, r)
    ωx, ωy = inv(σx^2), inv(σy^2)
    α = sqrt(ωx*ωy)
    Z = ωx*ωy / (b^2*ωy + ωx - 2*b*r*α)
    U, V = x - x̄, y - ȳ
    β = Z * (U/ωy + b*V/ωx - (b*U+V)*r/α)
    return Z, x̄ + β
end

# Correlation and covariance of two `Measurement`s
function paircov(a::Measurement{T}, b::Measurement{T}) where T
    overlap = keys(a.der) ∩ keys(b.der)
    return isempty(overlap) ? zero(T) : sum(overlap) do var
        a.der[var] * b.der[var] * var[2]^2
    end
end
paircor(a::Measurement, b::Measurement) = paircov(a,b)/(a.err*b.err)