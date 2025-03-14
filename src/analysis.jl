# Our overarching analysis type.
# Must contain a vector of means μ, standard deviations σ, and a covariance matrix Σ
abstract type Analysis{T<:AbstractFloat} <: Data{T} end

# Generic concrete implementation
struct BivariateAnalysis{T} <: Analysis{T}
    μ::SVector{2,T}
    σ::SVector{2,T}
    Σ::SMatrix{2,2,T,4}
end
function BivariateAnalysis(r1::Number, σ1::Number, r2::Number, σ2::Number, correlation::Number; T=Float64)
    cov = σ1 * σ2 * correlation
    Σ = SMatrix{2,2,T}(σ1^2, cov, cov, σ2^2)
    σ = SVector{2,T}(σ1, σ2)
    μ = SVector{2,T}(r1, r2)
    BivariateAnalysis(μ, σ, Σ)
end
BivariateAnalysis(μ::AbstractVector{T}, σ::AbstractVector, Σ::AbstractMatrix) where {T} = BivariateAnalysis{T}(SVector{2,T}(μ), SVector{2,T}(σ), SMatrix{2,2,T,4}(Σ))
BivariateAnalysis(μ::AbstractVector{T}, Σ::AbstractMatrix) where {T} = BivariateAnalysis{T}(SVector{2,T}(μ), SVector{2,T}(sqrt.(diag(Σ))), SMatrix{2,2,T,4}(Σ))
function BivariateAnalysis(x1::AbstractVector, x2::AbstractVector)
    μ1, μ2 = nanmean(x1), nanmean(x2)
    σ1, σ2 = nanstd(x1, mean=μ1), nanstd(x2, mean=μ2)
    σ12 = nancov(x1,x2)
    Σ = SMatrix{2,2}(σ1, σ12, σ12, σ2)
    BivariateAnalysis(SVector(μ1, μ2), SVector(σ1, σ2), Σ)
end


age(r::Number, λ::Number) = log(1+r)/λ
ratio(t::Number, λ::Number) = exp(λ*t) - 1

# Extend Base.isnan to return true if any component of the Analysis is NaN
Base.isnan(a::Analysis) = any(isnan, a.μ) || any(isnan, a.σ) || any(isnan, a.Σ)

# A moment in time
struct Age{T<:AbstractFloat} <: Data{T}
    mean::T
    sigma::T
end
Age(μ, σ) = Age(float(μ), float(σ))
Age(x) = Age(value(x), stdev(x))
Base.isless(x::Age, y::Age) = isless(x.mean, y.mean)

# A duration of time
struct Interval{T<:AbstractFloat} <: Data{T}
    min::T
    min_sigma::T
    max::T
    max_sigma::T
end
Interval(lμ, lσ, uμ, uσ) = Interval(float(lμ), float(lσ), float(uμ), float(uσ))
Interval(l, u) = Interval(value(l), stdev(l), value(u), stdev(u))
Base.min(x::Interval{T}) where {T} = Age{T}(x.min, x.min_sigma) 
Base.max(x::Interval{T}) where {T} = Age{T}(x.max, x.max_sigma)

# A confidence or credible interval with 95% bounds
struct CI{T<:AbstractFloat} <: Data{T}
    mean::T
    sigma::T
    median::T
    lower::T
    upper::T
end
function CI(x::AbstractVector{T}) where {T}
    xₜ = copy(x)
    Tₒ = float(T)
    mean = nanmean(xₜ)
    CI{Tₒ}(mean,
        nanstd(xₜ; mean),
        nanmedian!(xₜ),
        nanpctile!(xₜ, 2.5),
        nanpctile!(xₜ, 97.5),
    )
end
Base.isless(x::CI, y::CI) = isless(x.mean, y.mean)

# A type to hold a 2d covariance ellipse for any pair of measurements
struct Ellipse{T} <: Data{T}
    x::Vector{T}
    y::Vector{T}
    x₀::T
    y₀::T
    σx₀::T
    σy₀::T
end

# Make an ellipse from a Analysis object
function Ellipse(d::Analysis;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer=50,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return Ellipse(d, a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function Ellipse(d::Analysis, a, b, θ; npoints::Integer=50)
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

function datalimits(ellipses::Array{<:Ellipse})
    xmin = minimum(minimum.(x.(ellipses)))
    xmax = maximum(maximum.(x.(ellipses)))
    ymin = minimum(minimum.(y.(ellipses)))
    ymax = maximum(maximum.(y.(ellipses)))

    return xmin, xmax, ymin, ymax
end

datalimits(analyses::Array{<:Analysis}) = datalimits(Ellipse.(analyses))

x(e::Ellipse) = e.x
y(e::Ellipse) = e.y


# Convenience methods for possibly obtaining values or uncertainties
# Generic fallback methods for things that don't have uncertainties
value(x) = x
stdev(x::T) where {T} = zero(T)
# Specialized methods for `CI`s
value(x::CI{T}) where {T} = x.mean::T
stdev(x::CI{T}) where {T} = x.sigma::T
# Specialized methods for `Age`s
value(x::Age{T}) where {T} = x.mean::T
stdev(x::Age{T}) where {T} = x.sigma::T
# Specialized methods for `Measurement`s
value(x::Measurement{T}) where {T} = x.val::T
stdev(x::Measurement{T}) where {T} = x.err::T


# Deprecations of old methods
function val(x)
    @warn "`Isoplot.val` is deprecated, use `value` instead"
    return value(x)
end
function err(x)
    @warn "`Isoplot.err` is deprecated, use `stdev` instead"
    return stdev(x)
end