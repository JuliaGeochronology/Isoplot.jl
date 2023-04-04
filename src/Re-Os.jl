# Selby et al. (2007) 10.1016/j.gca.2007.01.008
λ187Re = 1.6689e-5 ± 0.0031e-5 # 1/Myr, calibrated against U-Pb. # 2-sigma?
export λ187Re

struct ReOsAnalysis{T} <: Analysis{T}
    μ::Vector{T}
    σ::Vector{T}
    Σ::Matrix{T}
end
ReOsAnalysis(μ::Vector{T}, σ::Vector{T}) where {T} = ReOsAnalysis{T}(μ, σ, diagm(σ).^2)
