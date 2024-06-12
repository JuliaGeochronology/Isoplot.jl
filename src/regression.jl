using SpecialFunctions: erfc

"""
Apply Chauvenet's criterion to a set of data to identify outliers.

The function calculates the z-scores of the data points, and then calculates the probability `p` of observing a value as extreme as the z-score under the assumption of normal distribution.
It then applies Chauvenet's criterion, marking any data point as an outlier if `2 * N * p < 1.0`, where `N` is the total number of data points.
"""
function chauvenet_func(μ::Vector{T}, σ::Vector) where {T}
    mean_val = mean(μ)
    N = length(μ)

    z_scores = abs.(μ .- mean_val) ./ σ
    p = 0.5 * erfc.(z_scores ./ sqrt(2.0))

    criterion = 2 * N * p
    selected_data = criterion .>= 1.0

    # add @info about number of outliers
    @info "Excluding $(N - sum(selected_data)) outliers based on Chauvenet's criterion."

    return selected_data
end

## --- Weighted means
"""
```julia
wμ, wσ, mswd = wmean(μ, σ; corrected=true, chauvenet=false)
wμ ± wσ, mswd = wmean(μ ± σ; corrected=true, chauvenet=false)
```
The weighted mean, with or without the "geochronologist's MSWD correction" to uncertainty.
You may specify your means and standard deviations either as separate vectors `μ` and `σ`,
or as a single vector `x` of `Measurement`s equivalent to `x = μ .± σ`

In all cases, `σ` is assumed to reported as _actual_ sigma (i.e., 1-sigma).

If `corrected=true`, the resulting uncertainty of the weighted mean is expanded by a factor
of `sqrt(mswd)` to attempt to account for dispersion dispersion when the MSWD is greater than `1`

If `chauvenet=true`, outliers will be removed before the computation of the weighted mean 
using Chauvenet's criterion.

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
(-0.32 ± 0.32, 0.8192171477885678)

julia> wmean(x .± y./10)
(-0.322 ± 0.032, 81.9217147788568)

julia> wmean(x .± y./10, corrected=true)
(-0.32 ± 0.29, 81.9217147788568)
```
"""
function wmean(μ::Collection1D{T}, σ::Collection1D{T}; corrected::Bool=true, chauvenet::Bool=false) where {T}
    if chauvenet
        not_outliers = chauvenet_func(μ, σ)
        μ = μ[not_outliers]
        σ = σ[not_outliers]
    end

    sum_of_values = sum_of_weights = χ² = zero(float(T))
    @inbounds for i in eachindex(μ,σ)
        σ² = σ[i]^2
        sum_of_values += μ[i] / σ²
        sum_of_weights += one(T) / σ²
    end
    wμ = sum_of_values / sum_of_weights

    @inbounds for i in eachindex(μ,σ)
        χ² += (μ[i] - wμ)^2 / σ[i]^2
    end
    mswd = χ² / (length(μ)-1)
    wσ = if corrected
        sqrt(max(mswd,1) / sum_of_weights)
    else
        sqrt(1 / sum_of_weights)
    end
    return wμ, wσ, mswd
end

function wmean(x::AbstractVector{Measurement{T}}; corrected::Bool=true, chauvenet::Bool=false) where {T}

    if chauvenet
        μ, σ = val.(x), err.(x)
        not_outliers = chauvenet_func(μ, σ)
        x = x[not_outliers]
    end

    wμ, wσ, mswd = wmean(val.(x), Measurements.cov(x); corrected)

    return wμ ± wσ, mswd
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
    corrected && (σ²ₓ̄ *= max(mswd,1))

    return x̄, sqrt(σ²ₓ̄), mswd
end

# Legacy methods, for backwards compatibility
awmean(args...) = wmean(args...; corrected=false)
gwmean(args...) = wmean(args...; corrected=true)


distwmean(x...; corrected::Bool=true) = distwmean(x; corrected)
function distwmean(x::NTuple{N, <:AbstractVector}; corrected::Bool=true) where {N}
    σₓ = vstd.(x)
    wₓ = 1 ./ σₓ.^2
    wₜ = sum(wₓ)
    c = similar(first(x))
    @inbounds for i in eachindex(x...)
        c[i] = sum(getindex.(x, i) .* wₓ)/wₜ
        if corrected
            c[i] += sqrt(sum(wₓ.*abs2.(getindex.(x,i).-c[i]))/(wₜ*(N-1))) * randn()
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
function mswd(μ::Collection{T}, σ::Collection; chauvenet=false) where {T}
    if chauvenet
        not_outliers = chauvenet_func(μ, σ)
        μ = μ[not_outliers]
        σ = σ[not_outliers]
    end

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

function mswd(x::AbstractVector{Measurement{T}}; chauvenet=false) where {T}

    if chauvenet
        not_outliers = chauvenet_func(val.(x), err.(x))
        x = x[not_outliers]
    end

    wμ, wσ, mswd = wmean(val.(x), Measurements.cov(x))

    return mswd
end

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
struct YorkFit{T<:Number}
    intercept::Measurement{T}
    slope::Measurement{T}
    xm::T
    ym::Measurement{T}
    mswd::T
end

"""
```julia
yorkfit(x, σx, y, σy, [r])
yorkfit(x::Vector{<:Measurement}, y::Vector{<:Measurement}, [r])
yorkfit(d::Vector{<:Analysis})
```
Uses the York (1968) two-dimensional least-squares fit to calculate `a`, `b`,
and uncertanties `σa`, `σb` for the equation `y = a + bx`, given `x`, `y`,
uncertaintes `σx`, and `σy`, and optially covarances `r`.

For further reference, see:
York, Derek (1968) "Least squares fitting of a straight line with correlated errors"
Earth and Planetary Science Letters 5, 320-324. doi: 10.1016/S0012-821X(68)80059-7

### Examples
```julia
julia> x = (1:100) .+ randn.();

julia> y = 2*(1:100) .+ randn.();

julia> yorkfit(x, ones(100), y, ones(100))
YorkFit{Float64}:
Least-squares linear fit of the form y = a + bx where
  intercept a : -0.29 ± 0.2 (1σ)
  slope b     : 2.0072 ± 0.0035 (1σ)
  MSWD        : 0.8136665223891004
```
"""
yorkfit(x::Vector{Measurement{T}}, y::Vector{Measurement{T}}, r=zero(T); iterations=10) where {T} = yorkfit(val.(x), err.(x), val.(y), err.(y), r; iterations)
function yorkfit(d::Collection{<:Analysis{T}}; iterations=10) where {T}
    # Using NTuples instead of Arrays here avoids allocations and should be
    # much more efficient for relatively small N, but could be less efficient
    # for large N (greater than ~100)
    x = ntuple(i->d[i].μ[1], length(d))
    y = ntuple(i->d[i].μ[2], length(d))
    σx = ntuple(i->d[i].σ[1], length(d))
    σy = ntuple(i->d[i].σ[2], length(d))
    r = ntuple(i->d[i].Σ[1,2], length(d))
    yorkfit(x, σx, y, σy, r; iterations)
end
function yorkfit(x, σx, y, σy, r=vcor(x,y); iterations=10)

    ## For an initial estimate of slope and intercept, calculate the
    # ordinary least-squares fit for the equation y=a+bx
    a, b = lsqfit(x, y)

    # Prepare for York fit
    ∅ = zero(float(eltype(x)))
    ωx = 1.0 ./ σx.^2           # x weights
    ωy = 1.0 ./ σy.^2           # y weights
    α = sqrt.(ωx .* ωy)

    ## Perform the York fit (must iterate)
    Z = @. ωx*ωy / (b^2*ωy + ωx - 2*b*r*α)

    x̄ = vsum(Z.*x) / vsum(Z)
    ȳ = vsum(Z.*y) / vsum(Z)

    U = x .- x̄
    V = y .- ȳ

    if Z isa NTuple
        Z = collect(Z)
        U = collect(U)
        V = collect(V)
    end

    sV = @. Z^2 * V * (U/ωy + b*V/ωx - r*V/α)
    sU = @. Z^2 * U * (U/ωy + b*V/ωx - b*r*U/α)
    b = vsum(sV) / vsum(sU)

    a = ȳ - b * x̄
    for _ in 2:iterations
        @. Z = ωx*ωy / (b^2*ωy + ωx - 2*b*r*α)

        ΣZ, ΣZx, ΣZy = ∅, ∅, ∅
        @inbounds for i in eachindex(Z)
            ΣZ += Z[i]
            ΣZx += Z[i] * x[i]
            ΣZy += Z[i] * y[i]
        end
        x̄ = ΣZx / ΣZ
        ȳ = ΣZy / ΣZ

        @. U = x - x̄
        @. V = y - ȳ

        @. sV = Z^2 * V * (U/ωy + b*V/ωx - r*V/α)
        @. sU = Z^2 * U * (U/ωy + b*V/ωx - b*r*U/α)
        b = sum(sV) / sum(sU)

        a = ȳ - b * x̄
    end

    ## 4. Calculate uncertainties and MSWD
    β = @. Z * (U/ωy + b*V/ωx - (b*U+V)*r/α)

    u = x̄ .+ β
    v = ȳ .+ b.*β

    xm = vsum(Z.*u)./vsum(Z)
    ym = vsum(Z.*v)./vsum(Z)

    σb = sqrt(1.0 ./ vsum(Z .* (u .- xm).^2))
    σa = sqrt(1.0 ./ vsum(Z) + xm.^2 .* σb.^2)
    σym = sqrt(1.0 ./ vsum(Z))

    # MSWD (reduced chi-squared) of the fit
    mswd = 1.0 ./ length(x) .* vsum(@. (y - a - b*x)^2 / (σy^2 + b^2 * σx^2) )

    ## Results
    return YorkFit(a ± σa, b ± σb, xm, ym ± σym, mswd)
end
