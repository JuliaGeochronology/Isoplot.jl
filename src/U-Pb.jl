# Decay constants:
const λ238U = log(2)/(4.4683e3 ± 0.0024e3) # Jaffey, 1/Myr
const λ235U = 9.8569E-4 ± 0.0017E-4 # Schoene, 1/Myr [U-238 years]
export λ238U, λ235U

const λ235U_jaffey = log(2)/(7.0381e2 ± 0.0048e2) # Jaffey, 1/Myr
export λ235U_jaffey

"""
```
struct UPbAnalysis{T} <: Analysis{T}
```
Core type for U-Pb analyses.
Has fields
```
μ :: Vector{T<:AbstractFloat}
σ :: Vector{T<:AbstractFloat}
Σ :: Matrix{T<:AbstractFloat}
```
where `μ` contains the means
```
μ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
```
where `σ` contains the standard deviations
```
σ = [σ²⁰⁷Pb²³⁵U, σ²⁰⁶Pb²³⁸U]
```
and Σ contains the covariance matrix
```
Σ = [σ₇_₅^2 σ₇_₅*σ₃_₈
     σ₇_₅*σ₃_₈ σ₃_₈^2]
```
If `σ` is not provided, it will be automatically calculated from `Σ`,
given that `σ.^2 = diag(Σ)`.
"""
struct UPbAnalysis{T} <: Analysis{T}
    μ::Vector{T}
    σ::Vector{T}
    Σ::Matrix{T}
end

"""
```julia
UPbAnalysis(r²⁰⁷Pb²³⁵U, σ²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U, σ²⁰⁶Pb²³⁸U, correlation; T=Float64)
```
Construct a `UPbAnalysis` object from individual isotope ratios and (1-sigma!) uncertainties.

### Examples
```
julia> UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183)
UPbAnalysis{Float64}([22.6602, 0.40864], [0.00030625000000000004 2.4746942500000003e-6; 2.4746942500000003e-6 2.8900000000000004e-8])
```
"""
function UPbAnalysis(r²⁰⁷Pb²³⁵U::Number, σ²⁰⁷Pb²³⁵U::Number, r²⁰⁶Pb²³⁸U::Number, σ²⁰⁶Pb²³⁸U::Number, correlation::Number; T=Float64)
    cov = σ²⁰⁷Pb²³⁵U * σ²⁰⁶Pb²³⁸U * correlation
    Σ = T[σ²⁰⁷Pb²³⁵U^2  cov
          cov  σ²⁰⁶Pb²³⁸U^2]
    σ = T[σ²⁰⁷Pb²³⁵U, σ²⁰⁶Pb²³⁸U]
    μ = T[r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
    UPbAnalysis(μ, σ, Σ)
end
UPbAnalysis(μ::Vector{T}, Σ::Matrix{T}) where {T} = UPbAnalysis{T}(μ, sqrt.(diag(Σ)), Σ)

# 75 and 68 ages
function age(d::UPbAnalysis)
    a75 = log(1 + d.μ[1] ± d.σ[1])/λ235U
    a68 = log(1 + d.μ[2] ± d.σ[2])/λ238U
    return a75, a68
end

# Percent discordance
function discordance(d::UPbAnalysis)
    μ75 = log(1 + d.μ[1])/λ235U.val
    μ68 = log(1 + d.μ[2])/λ238U.val
    return (μ75 - μ68) / μ75 * 100
end

# Add custom methods to Base.rand to sample from a UPbAnalysis
Base.rand(d::UPbAnalysis) = rand(MvNormal(d.μ, d.Σ))
Base.rand(d::UPbAnalysis, n::Integer) = rand(MvNormal(d.μ, d.Σ), n)
Base.rand(d::UPbAnalysis, dims::Dims) = rand(MvNormal(d.μ, d.Σ), dims)
