λ147Sm = 6.524-6 ± 0.024e-6/2 # 1/Myr, Villa et al. (2020) 10.1016/j.gca.2020.06.022
λ146Sm = log(2)/(103.1 ± 4.5) # 1/Myr, Meissner et al. (1987)
export λ147Sm, λ146Sm

struct SmNdAnalysis{T} <: ParentDaughterAnalysis{T}
    data::Analysis2D{T}
end
SmNdAnalysis(args...; kwargs...) = SmNdAnalysis(Analysis(args...; kwargs...))

lambda(::Type{<:SmNdAnalysis}) = λ147Sm