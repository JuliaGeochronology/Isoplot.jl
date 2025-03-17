# Nebel et al. (2011) 10.1016/j.epsl.2010.11.004
λ87Rb = 1.393e-5 ± 0.004e-5/2 # 1/Myr, calibrated against U-Pb
export λ87Rb

struct RbSrAnalysis{T} <: ParentDaughterAnalysis{T}
    data::Analysis2D{T}
end
RbSrAnalysis(args...; kwargs...) = RbSrAnalysis(Analysis(args...; kwargs...))
