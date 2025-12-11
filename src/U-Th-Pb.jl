# Decay constants:
const λ238U = log(2)/(4.4683e3 ± 0.0024e3) # [1/Myr] Jaffey et al. (1971)
const λ235U = 9.8569E-4 ± 0.0110E-4/2 # [1/Myr] Schoene et al. (2006)
const λ235U_internal = 9.8569E-4 ± 0.0017E-4/2 # [1/Myr] Schoene et al. (2006), including only internal uncertainty [U-238 years]
export λ238U, λ235U

const λ235U_jaffey = log(2)/(7.0381e2 ± 0.0048e2) # [1/Myr] Jaffey et al (1971)
export λ235U_jaffey

const λ232Th = log(2)/(1.401e4 ± 0.007e4) # [1/Myr] LeRoux & Glendenin (1963)
export λ232Th


"""
```
struct UThPbAnalysis{T} <: NoInitialDaughterAnalysis{T}
```
Core type for U-Th-Pb analyses.
Wraps an Analysis object which has fields
```
μ :: SVector{T<:AbstractFloat}
σ :: SVector{T<:AbstractFloat}
Σ :: SMatrix{T<:AbstractFloat}
```
where `μ` contains the means
```
μ = [μ²⁰⁷Pb/²³⁵U, μ²⁰⁶Pb/²³⁸U, μ²⁰⁸Pb/²³²Th]
```
where `σ` contains the standard deviations
```
σ = [σ²⁰⁷Pb²³⁵U, σ²⁰⁶Pb²³⁸U, σ²⁰⁸Pb/²³²Th]
```
and Σ contains the covariance matrix
```
Σ = [ σ₇_₅^2  σ₇_₅σ₃_₈ σ₇_₅σ₈_₂
     σ₇_₅σ₃_₈  σ₃_₈^2  σ₃_₈σ₈_₂
     σ₇_₅σ₈_₂ σ₃_₈σ₈_₂  σ₈_₂^2 ]
```
If `σ` is not provided, it will be automatically calculated from `Σ`,
given that `σ.^2 = diag(Σ)`.
"""
struct UThPbAnalysis{T} <: NoInitialDaughterAnalysis{T}
    data::Analysis3D{T}
end


"""
```julia
UThPbAnalysis(μ, Σ)
```
Construct a `UPbAnalysis` object from vector  of isotope ratios and (1-sigma!) covariance matrix.
"""
UThPbAnalysis(args...; kwargs...) = UThPbAnalysis(Analysis(args...; kwargs...))


"""
```
struct UPbAnalysis{T} <: NoInitialDaughterAnalysis{T}
```
Core type for U-Pb analyses.
Wraps an Analysis object which has fields
```
μ :: SVector{T<:AbstractFloat}
σ :: SVector{T<:AbstractFloat}
Σ :: SMatrix{T<:AbstractFloat}
```
where `μ` contains the means
```
μ = [μ²⁰⁷Pb/²³⁵U, μ²⁰⁶Pb/²³⁸U]
```
where `σ` contains the standard deviations
```
σ = [σ²⁰⁷Pb/²³⁵U, σ²⁰⁶Pb/²³⁸U]
```
and Σ contains the covariance matrix
```
Σ = [ σ₇_₅^2  σ₇_₅σ₃_₈
     σ₇_₅σ₃_₈  σ₃_₈^2 ]
```
If `σ` is not provided, it will be automatically calculated from `Σ`,
given that `σ.^2 = diag(Σ)`.
"""
struct UPbAnalysis{T} <: NoInitialDaughterAnalysis{T}
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
UPbAnalysis(args...; kwargs...) = UPbAnalysis(Analysis(args...; kwargs...))


# Pb-207/U-235 and Pb-206/U-238, and Pb-208/Th-232 ages
function age(d::UPbAnalysis)
    return age75(d), age68(d)
end
function age(d::UThPbAnalysis)
    return age75(d), age68(d), age82(d)
end

function age68(d::Union{UPbAnalysis{T}, UThPbAnalysis{T}}) where {T}
    return log(1 + max(mean(d)[2], zero(T)) ± std(d)[2])/λ238U
end
function age75(d::Union{UPbAnalysis{T}, UThPbAnalysis{T}}) where {T}
    return log(1 + max(mean(d)[1], zero(T)) ± std(d)[1])/λ235U
end
function age82(d::UThPbAnalysis{T}) where {T}
    return log(1 + max(mean(d)[3], zero(T)) ± std(d)[3])/λ232Th
end

# Percent discordance
function discordance(d::UPbAnalysis{T}) where {T}
    μ75 = mean(d)[1] > 0  ? log(1 + mean(d)[1])/value(λ235U) : T(NaN)
    μ68 = mean(d)[2] > 0  ? log(1 + mean(d)[2])/value(λ238U) : T(NaN)
    return (μ75 - μ68) / μ75 * 100
end

function stacey_kramers(t)
    if 3700 <= t < 4570
        t0 = 3700
        r64 = 11.152
        r74 = 12.998
        r84 = 31.230
        U_Pb = 7.19
        Th_Pb = 32.21
    elseif t < 3700
        t0 = 0
        r64 = 18.700
        r74 = 15.628
        r84 = 38.630
        U_Pb = 9.74
        Th_Pb = 37.19
    else
        t0 = NaN
        r64 = NaN
        r74 = NaN
        r84 = NaN
        U_Pb = NaN
        Th_Pb = NaN
    end

    r64 -= ((exp(value(λ238U)*t)-1) - (exp(value(λ238U)*t0)-1)) * U_Pb
    r74 -= ((exp(value(λ238U)*t)-1) - (exp(value(λ238U)*t0)-1)) * U_Pb/137.818
    r84 -= ((exp(value(λ232Th)*t)-1) - (exp(value(λ232Th)*t0)-1)) * Th_Pb

    return r64, r74, r84
end
