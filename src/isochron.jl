# Converting from age to ratio and back
age(r::Number, λ::Number) = log(1+r)/λ
ratio(t::Number, λ::Number) = exp(λ*t) - 1

# Define isochron type
struct Isochron{T,A<:ParentDaughterAnalysis{T}} <: Data{T}
    line::YorkFit{T}
end

# Calculate an isochron
isochron(x::Collection{A}) where {T,A<:ParentDaughterAnalysis{T}} = Isochron{T,A}(yorkfit(x))

# Calculate the age of an isochron
age(x::Isochron{T,A}, λ=lambda(A)) where {T,A} = age(x.line.slope, λ)



