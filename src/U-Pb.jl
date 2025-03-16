# Decay constants:
const λ238U = log(2)/(4.4683e3 ± 0.0024e3) # Jaffey, 1/Myr
const λ235U = 9.8569E-4 ± 0.0110E-4/2 # Schoene, 1/Myr
const λ235U_internal = 9.8569E-4 ± 0.0017E-4/2 # Schoene, 1/Myr, including only internal uncertainty [U-238 years]
export λ238U, λ235U

const λ235U_jaffey = log(2)/(7.0381e2 ± 0.0048e2) # Jaffey, 1/Myr
export λ235U_jaffey

"""
```
struct UPbAnalysis{T} <: Analysis{T}
```
Core type for U-Pb analyses.
Has fields
```
μ :: Vector{T<:AbstractFloat}
σ :: Vector{T<:AbstractFloat}
Σ :: Matrix{T<:AbstractFloat}
```
where `μ` contains the means
```
μ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]
```
where `σ` contains the standard deviations
```
σ = [σ²⁰⁷Pb²³⁵U, σ²⁰⁶Pb²³⁸U]
```
and Σ contains the covariance matrix
```
Σ = [σ₇_₅^2 σ₇_₅*σ₃_₈
     σ₇_₅*σ₃_₈ σ₃_₈^2]
```
If `σ` is not provided, it will be automatically calculated from `Σ`,
given that `σ.^2 = diag(Σ)`.
"""
struct UPbAnalysis{T} <: AbstractAnalysis{T}
    data::Analysis2D{T}
end


"""
```julia
UPbAnalysis(r²⁰⁷Pb²³⁵U, σ²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U, σ²⁰⁶Pb²³⁸U, correlation; T=Float64)
```
Construct a `UPbAnalysis` object from individual isotope ratios and (1-sigma!) uncertainties.

### Examples
```
julia> UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183)
UPbAnalysis{Float64}([22.6602, 0.40864], [0.00030625000000000004 2.4746942500000003e-6; 2.4746942500000003e-6 2.8900000000000004e-8])
```
"""
UPbAnalysis(args...) = UPbAnalysis(Analysis(args...))


# 75 and 68 ages
function age(d::UPbAnalysis)
    a75 = log(1 + mean(d)[1] ± std(d)[1])/λ235U
    a68 = log(1 + mean(d)[2] ± std(d)[2])/λ238U
    return a75, a68
end

function age68(d::UPbAnalysis)
    log(1 + mean(d)[2] ± std(d)[2])/λ238U
end
function age75(d::UPbAnalysis)
    log(1 + mean(d)[1] ± std(d)[1])/λ235U
end
# Percent discordance
function discordance(d::UPbAnalysis)
    μ75 = log(1 + mean(d)[1])/value(λ235U)
    μ68 = log(1 + mean(d)[2])/value(λ238U)
    return (μ75 - μ68) / μ75 * 100
end

# Add custom methods to Base.rand to sample from a UPbAnalysis
Base.rand(d::UPbAnalysis) = rand(MvNormal(mean(d), cov(d)))
Base.rand(d::UPbAnalysis, n::Integer) = rand(MvNormal(mean(d), cov(d)), n)
Base.rand(d::UPbAnalysis, dims::Dims) = rand(MvNormal(mean(d), cov(d)), dims)

function stacey_kramers(t)
    if 3700 <= t < 4570
        t0 = 3700
        r64 = 11.152
        r74 = 12.998
        U_Pb = 7.19
    elseif t < 3700
        t0 = 0
        r64 = 18.700
        r74 = 15.628
        U_Pb = 9.74
    else
        t0 = NaN
        r64 = NaN
        r74 = NaN
        U_Pb = NaN
    end

    r64 -= ((exp(value(λ238U)*t)-1) - (exp(value(λ238U)*t0)-1)) * U_Pb
    r74 -= ((exp(value(λ238U)*t)-1) - (exp(value(λ238U)*t0)-1)) * U_Pb/137.818

    return r64, r74
end
