λ147Sm = 6.524-6 ± 0.024e-6/2 # 1/Myr, Villa et al. (2020) 10.1016/j.gca.2020.06.022
λ146Sm = log(2)/(103.1 ± 4.5) # 1/Myr, Meissner et al. (1987)
export λ147Sm, λ147Sm

struct SmNdAnalysis{T} <: Analysis{T}
    μ::Vector{T}
    σ::Vector{T}
    Σ::Matrix{T}
end
SmNdAnalysis(μ::Vector{T}, σ::Vector{T}) where {T} = SmNdAnalysis{T}(μ, σ, diagm(σ).^2)
