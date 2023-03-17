# Decay constants:
# const λ235U = log(2)/(7.0381*10^2) ± log(2)/(7.0381*10^2)*0.0048/7.0381 # Jaffey, 1/Myr
const λ238U = log(2)/(4.4683*10^3) ± log(2)/(4.4683*10^3)*0.0024/4.4683 # Jaffey, 1/Myr
const λ235U = 9.8569E-4 ± 0.0017E-4 # Schoene, 1/Myr [U-238 years]
export λ238U, λ235U

"""
```
struct UPbAnalysis{T<:AbstractFloat}
```
Core type for U-Pb analyses.
Has fields
```
μ :: Vector{T<:AbstractFloat}
Σ :: Matrix{T<:AbstractFloat}
```
where `μ` contains the means
```
μ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
```
and Σ contains the covariance matrix
```
Σ = [σ₇_₅^2 σ₇_₅*σ₃_₈
     σ₇_₅*σ₃_₈ σ₃_₈^2]
```
"""
struct UPbAnalysis{T<:AbstractFloat}
    μ::Vector{T}
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
    μ = T[r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
    UPbAnalysis(μ, Σ)
end


# Add custom methods to Base.rand to sample from a UPbAnalysis
Base.rand(d::UPbAnalysis) = rand(MvNormal(d.μ, d.Σ))
Base.rand(d::UPbAnalysis, n::Integer) = rand(MvNormal(d.μ, d.Σ), n)
Base.rand(d::UPbAnalysis, dims::Dims) = rand(MvNormal(d.μ, d.Σ), dims)
