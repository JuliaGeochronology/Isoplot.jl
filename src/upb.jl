# Decay constants
# Uranium 238 decay constant and uncertainty
const λ238U = log(2)/(4.4683*10^9) ± log(2)/(4.4683*10^9)*0.0024/4.4683 # Jaffey, 1/Years
export λ238U

# Uranium 235 decay constant and uncertainty
# const λ235U = log(2)/(7.0381*10^8) ± log(2)/(7.0381*10^8)*0.0048/7.0381 # Jaffey, 1/Years
const λ235U = 9.8569E-10 ± 0.0017E-10 # Schoene, 1/Years [U238 years]
export λ235U

# Core type and constructors for U-Pb analyses
# contains
# μ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
# Σ = [σ₇_₅^2 σ₇_₅*σ₃_₈
#      σ₇_₅*σ₃_₈ σ₃_₈^2]
#
struct UPbAnalysis{T<:AbstractFloat}
    μ::Vector{T}
    Σ::Matrix{T}
end
function UPbAnalysis(r²⁰⁷Pb²³⁵U::Number, σ²⁰⁷Pb²³⁵U::Number, r²⁰⁶Pb²³⁸U::Number, σ²⁰⁶Pb²³⁸U::Number, correlation::Number; T=Float64)
    cov = σ²⁰⁷Pb²³⁵U * σ²⁰⁶Pb²³⁸U * correlation
    Σ = T[σ²⁰⁷Pb²³⁵U^2  cov
          cov  σ²⁰⁶Pb²³⁸U^2]
    μ = T[r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
    UPbAnalysis(μ, Σ)
end
export UPbAnalysis

# Add custom methods to Base.rand to sample from a UPbAnalysis
Base.rand(d::UPbAnalysis) = rand(MvNormal(d.μ, d.Σ))
Base.rand(d::UPbAnalysis, n::Integer) = rand(MvNormal(d.μ, d.Σ), n)
Base.rand(d::UPbAnalysis, dims::Dims) = rand(MvNormal(d.μ, d.Σ), dims)
