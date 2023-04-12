## --- Weighted means

"""
```julia
wμ, wσ, mswd = wmean(μ, σ; corrected=false)
wμ ± wσ, mswd = wmean(μ ± σ; corrected=false)
```
The weighted mean, with or without the "geochronologist's MSWD correction" to uncertainty.
You may specify your means and standard deviations either as separate vectors `μ` and `σ`,
or as a single vector `x` of `Measurement`s equivalent to `x = μ .± σ`

In all cases, `σ` is assumed to reported as _actual_ sigma (i.e., 1-sigma).

If `corrected=true`, the resulting uncertainty of the weighted mean is corrected
for dispersion when the MSWD is greater than `1` by multiplying by the square
root of the MSWD.

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
function wmean(μ::Collection{T}, σ::Collection; corrected::Bool=true) where {T}
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
function wmean(x::Collection{Measurement{T}}; corrected::Bool=true) where {T}
    sum_of_values = sum_of_weights = χ² = zero(float(T))
    @inbounds for i in eachindex(x)
        μ, σ² = val(x[i]), err(x[i])^2
        sum_of_values += μ / σ²
        sum_of_weights += one(T) / σ²
    end
    wμ = sum_of_values / sum_of_weights

    @inbounds for i in eachindex(x)
        μ, σ = val(x[i]), err(x[i])
        χ² += (μ - wμ)^2 / σ^2
    end
    mswd = χ² / (length(x)-1)
    wσ = if corrected
        sqrt(max(mswd,1) / sum_of_weights)
    else
        sqrt(1 / sum_of_weights)
    end
    return wμ ± wσ, mswd
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
function mswd(μ::Collection{T}, σ::Collection) where {T}
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

function mswd(x::Collection{Measurement{T}}) where {T}
    sum_of_values = sum_of_weights = χ² = zero(float(T))

    @inbounds for i in eachindex(x)
        w = 1 / err(x[i])^2
        sum_of_values += w * val(x[i])
        sum_of_weights += w
    end
    wx = sum_of_values / sum_of_weights

    @inbounds for i in eachindex(x)
        χ² += (val(x[i]) - wx)^2 / err(x[i])^2
    end

    return χ² / (length(x)-1)
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
function lsqfit(x::AbstractVector{T}, y::AbstractVector{<:Number}) where {T<:Number}
    A = similar(x, length(x), 2)
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
function yorkfit(d::Vector{<:Analysis{T}}; iterations=10) where {T}
    x, y, σx, σy, r = similar(d, T), similar(d, T), similar(d, T), similar(d, T), similar(d, T)
    @inbounds for i in eachindex(d)
        dᵢ = d[i]
        x[i], y[i] = dᵢ.μ[1], dᵢ.μ[2]
        σx[i], σy[i] = dᵢ.σ[1], dᵢ.σ[2]
        r[i] = dᵢ.Σ[2,1]
    end
    yorkfit(x, σx, y, σy, r; iterations)
end
function yorkfit(x, σx, y, σy, r=vcor(x,y); iterations=10)

    ## Check for and exclude missing data
    t = (x.==x) .& (y.==y) .& (σx.==σx) .& (σy.==σy) .& (r.==r)
    if any(t)
        x, y = x[t], y[t]
        σx, σy = σx[t], σy[t]
        r isa Vector && (r = r[t])
    end

    ## For an initial estimate of slope and intercept, calculate the
    # ordinary least-squares fit for the equation y=a+bx
    a, b = lsqfit(x, y)

    # Prepare for York fit
    ∅ = zero(float(eltype(x)))
    ωx = 1.0 ./ σx.^2           # x weights
    ωy = 1.0 ./ σy.^2           # y weights
    α = sqrt.(ωx .* ωy)

    ## Perform the York fit (must iterate)
    W = ωx.*ωy ./ (b^2*ωy + ωx - 2*b*r.*α)

    X̄ = vsum(W.*x) / vsum(W)
    Ȳ = vsum(W.*y) / vsum(W)

    U = x .- X̄
    V = y .- Ȳ

    sV = @. W^2 * V * (U/ωy + b*V/ωx - r*V/α)
    sU = @. W^2 * U * (U/ωy + b*V/ωx - b*r*U/α)
    b = vsum(sV) / vsum(sU)

    a = @. Ȳ - b * X̄
    for _ in 2:iterations
        @. W = ωx*ωy / (b^2*ωy + ωx - 2*b*r*α)

        ΣW, ΣWx, ΣWy = ∅, ∅, ∅
        @inbounds for i in eachindex(W,x,y)
            ΣW += W[i]
            ΣWx += W[i] * x[i]
            ΣWy += W[i] * y[i]
        end
        X̄ = ΣWx / ΣW
        Ȳ = ΣWy / ΣW

        @. U = x - X̄
        @. V = y - Ȳ

        @. sV = W^2 * V * (U/ωy + b*V/ωx - r*V/α)
        @. sU = W^2 * U * (U/ωy + b*V/ωx - b*r*U/α)
        b = sum(sV) / sum(sU)

        a = Ȳ - b .* X̄
    end

    ## 4. Calculate uncertainties and MSWD
    β = @. W * (U/ωy + b*V/ωx - (b*U+V)*r/α)

    u = X̄ .+ β
    v = Ȳ .+ b.*β

    xm = vsum(W.*u)./vsum(W)
    ym = vsum(W.*v)./vsum(W)

    σb = sqrt(1.0 ./ vsum(W .* (u .- xm).^2))
    σa = sqrt(1.0 ./ vsum(W) + xm.^2 .* σb.^2)

    # MSWD (reduced chi-squared) of the fit
    mswd = 1.0 ./ length(x) .* vsum( (y .- a.-b.* x).^2 ./ (σy.^2 + b.^2 .* σx.^2) )

    ## Results
    return YorkFit(a ± σa, b ± σb, mswd)
end
