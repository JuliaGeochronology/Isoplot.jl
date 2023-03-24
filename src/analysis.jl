# Our overarching analysis type.
# Must contain a vector of means μ, standard deviations σ, and a covariance matrix Σ
abstract type Analysis{T<:Float64} end

age(r::Number, λ::Number) = log(1+r)/λ
ratio(t::Number, λ::Number) = exp(λ*t) - 1

# Extend Base.isnan to return true if any component of the Analysis is NaN
Base.isnan(a::Analysis) = any(isnan, a.μ) || any(isnan, a.σ) || any(isnan, a.Σ)

# A type to hold a 2d covariance ellipse for any pair of measurements
struct Ellipse{T}
    x::Vector{T}
    y::Vector{T}
    x₀::T
    y₀::T
    σx₀::T
    σy₀::T
end

# Make an ellipse from a Analysis object
function ellipse(d::Analysis;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer=50,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return ellipse(d, a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function ellipse(d::Analysis, a, b, θ; npoints::Integer=50)
    x₀, y₀ = d.μ[1], d.μ[2]
    t = range(0, 2π, length=npoints)
    x = a*cos(θ)*cos.(t) .- b*sin(θ)*sin.(t) .+ x₀
    y = a*sin(θ)*cos.(t) .+ b*cos(θ)*sin.(t) .+ y₀
    return Ellipse(x, y, x₀, y₀, d.σ[1], d.σ[2])
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

# Convenience methods for possibly obtaining values or uncertainties
# Generic fallback methods for things that don't have uncertainties
val(x) = x
err(x::T) where {T} = zero(T)
# Specialized methods for `Measurement`s
val(m::Measurement{T}) where {T} = m.val::T
err(m::Measurement{T}) where {T} = m.err::T
