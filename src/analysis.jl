# Our overarching analysis type.
abstract type AbstractAnalysis{T<:AbstractFloat} <: Data{T} end
# Instances other than `Analysis` are intended to be wrapper types, and
# must contain an `Analysis` object in a field called `data`

# Concrete implementation
struct Analysis{N,T,N2} <: AbstractAnalysis{T}
    μ::SVector{N,T}
    σ::SVector{N,T}
    Σ::SMatrix{N,N,T,N2}
end
const Analysis1D{T} = Analysis{1,T,1}
const Analysis2D{T} = Analysis{2,T,4}
const Analysis3D{T} = Analysis{3,T,9}
const Analysis4D{T} = Analysis{4,T,16}
const Analysis5D{T} = Analysis{5,T,25}


function Analysis(x1::AbstractVector, x2::AbstractVector)
    @assert eachindex(x1) == eachindex(x2) "Vectors x1 and x2 must be of equal dimensions"
    μ1, μ2 = nanmean(x1), nanmean(x2)
    σ1, σ2 = nanstd(x1, mean=μ1), nanstd(x2, mean=μ2)
    σ12 = nancov(x1,x2)
    Σ = SMatrix{2,2}(σ1, σ12, σ12, σ2)
    Analysis(SVector(μ1, μ2), SVector(σ1, σ2), Σ)
end
function Analysis(r1::Number, σ1::Number, r2::Number, σ2::Number, correlation::Number; T=Float64)
    cov = σ1 * σ2 * correlation
    Σ = SMatrix{2,2,T}(σ1^2, cov, cov, σ2^2)
    σ = SVector{2,T}(σ1, σ2)
    μ = SVector{2,T}(r1, r2)
    return Analysis2D{T}(μ, σ, Σ)
end
function Analysis(μ::AbstractVector{T}, σ::AbstractVector, Σ::AbstractMatrix) where {T}
    N, N2 = length(μ), length(μ)*length(μ)
    return Analysis{N,T,N2}(SVector{N,T}(μ), SVector{N,T}(σ), SMatrix{N,N,T,N2}(Σ))
end
Analysis(μ::AbstractVector, Σ::AbstractMatrix) = Analysis(μ, sqrt.(diag(Σ)), Σ)

# Extend Base.isnan to return true if any component of the Analysis is NaN
Base.isnan(a::Analysis) = any(isnan, a.μ) || any(isnan, a.σ) || any(isnan, a.Σ)
Base.isnan(a::AbstractAnalysis) = isnan(a.data)

# Extend mean, std, and cov for Analysis objects
Distributions.mean(a::Analysis) = a.μ
Distributions.mean(a::AbstractAnalysis) = mean(a.data)
Distributions.std(a::Analysis) = a.σ
Distributions.std(a::AbstractAnalysis) = std(a.data)
Distributions.cov(a::Analysis) = a.Σ
Distributions.cov(a::AbstractAnalysis) = cov(a.data)


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
Ellipse(d::AbstractAnalysis, args...; kwargs...) = Ellipse(d.data, args...; kwargs...)
function Ellipse(d::Analysis2D;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer=50,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return Ellipse(d, a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function Ellipse(d::Analysis2D, a, b, θ; npoints::Integer=50)
    x₀, y₀ = d.μ[1], d.μ[2]
    t = range(0, 2π, length=npoints)
    x = a*cos(θ)*cos.(t) .- b*sin(θ)*sin.(t) .+ x₀
    y = a*sin(θ)*cos.(t) .+ b*cos(θ)*sin.(t) .+ y₀
    return Ellipse(x, y, x₀, y₀, std(d)[1], std(d)[2])
end

# Non-exported function: return semimajor and minor axes for a given U-Pb analysis
function ellipseparameters(d::Analysis2D{T}, sigmalevel::Number) where T

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

xvals(e::Ellipse) = e.x
yvals(e::Ellipse) = e.y

function datalimits(ellipses::Array{<:Ellipse})
    xmin = minimum(minimum.(xvals.(ellipses)))
    xmax = maximum(maximum.(xvals.(ellipses)))
    ymin = minimum(minimum.(yvals.(ellipses)))
    ymax = maximum(maximum.(yvals.(ellipses)))

    return xmin, xmax, ymin, ymax
end

datalimits(analyses::Array{<:AbstractAnalysis}) = datalimits(Ellipse.(analyses))


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