# Our overarching analysis type.
# Must contain a vector of means μ and a covariance matrix Σ
abstract type Analysis{T<:Float64} end

struct Ellipse{T}
    x::Vector{T}
    y::Vector{T}
    x₀::T
    y₀::T
end

# Make an ellipse from a Analysis object
function ellipse(d::Analysis;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer=50,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return ellipse(d.μ[1], d.μ[2], a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function ellipse(x₀, y₀, a, b, θ; npoints::Integer=50)
    t = range(0, 2π, length=npoints)
    x = a*cos(θ)*cos.(t) .- b*sin(θ)*sin.(t) .+ x₀
    y = a*sin(θ)*cos.(t) .+ b*cos(θ)*sin.(t) .+ y₀
    return Ellipse(x, y, x₀, y₀,)
end

# Non-exported function: return semimajor and minor axes for a given U-Pb analysis
function ellipseparameters(d::Analysis{T}, sigmalevel::Number) where T

    # Quickly exit if any NaNs
    any(isnan, d.Σ) && return T.((NaN, NaN, NaN))

    # Calculate eigenvectors and eigenvalues from the covariance matrix.
    # V: matrix of eigenvectors, D: diagonal matrix of eigenvalues
    F = eigen(d.Σ)
    # Find index of semimajor and semiminor axes
    major = argmax(F.values)
    minor = argmin(F.values)
    v = view(F.vectors, :, major)

    # Calculate angle of major axis of ellipse from horizontal
    θ = atan(v[2]/v[1])

    # Calculate length of semimajor and semiminor axes for given p-value
    a = T(sigmalevel)*sqrt(abs(F.values[major]))
    b = T(sigmalevel)*sqrt(abs(F.values[minor]))

    return a, b, θ
end

# Convenience methods for converting between Measurements and values and errors
@inline val(m::Measurement{T}) where {T} = m.val::T
@inline err(m::Measurement{T}) where {T} = m.err::T
val(x::Vector{<:Measurement}) = [val(m) for m in x]
err(x::Vector{<:Measurement}) = [err(m) for m in x]
