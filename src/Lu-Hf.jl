# Söderlund et al. (2004) 0.1016/S0012-821X(04)00012-3
λ176Lu = 1.867e-5 ± 0.008e-5 # 1/Myr, calibrated against U-Pb # 2-sigma?
export λ176Lu

struct LuHfAnalysis{T} <: Analysis{T}
    μ::Vector{T}
    σ::Vector{T}
    Σ::Matrix{T}
end
LuHfAnalysis(μ::Vector{T}, σ::Vector{T}) where {T} = LuHfAnalysis{T}(μ, σ, diagm(σ).^2)
