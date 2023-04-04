# Cheng et al. (2013) 10.1016/j.epsl.2013.04.006
const λ234U = log(2)/(245620e-6 ± 260e-6/2) # 1/Myr
const λ230Th = log(2)/(75584e-6 ± 110e-6/2) # 1/Myr
export λ234U, λ230Th

struct UThAnalysis{T} <: Analysis{T}
    μ::Vector{T}
    σ::Vector{T}
    Σ::Matrix{T}
end
UThAnalysis(μ::Vector{T}, σ::Vector{T}) where {T} = UThAnalysis{T}(μ, σ, diagm(σ).^2)
