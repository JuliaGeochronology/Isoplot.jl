module Isoplot

    using Reexport
    @reexport using NaNStatistics
    using LoopVectorization: @turbo
    using StatGeochemBase: inpolygon
    using DelimitedFiles: readdlm
    using StaticArrays: SVector, SMatrix, @SVector
    using Rotations: RotMatrix
    using LogExpFunctions
    using LinearAlgebra
    using Distributions
    using Measurements
    
    # Some utilities to support use of NTuples within functions that are expecting arrays
    const Collection{T} = Union{AbstractArray{T}, NTuple{N,T}} where N
    const Collection1D{T} = Union{AbstractVector{T}, NTuple{N,T}} where N
    arraylike(args...) = similar(args...)
    arraylike(x::NTuple{N,T}) where {N,T} = similar(Array{T}, N)
    arraylike(x::NTuple{N,T}, args...) where {N,T} = similar(Array{T}, args...)
    arraylike(x::NTuple{N}, ::Type{T}) where {N,T} = similar(Array{T}, N)
    arraylike(x::NTuple, ::Type{T}, args...) where {T} = similar(Array{T}, args...)

    # Age of Earth and the Solar System
    const t🜨 = 4.567e3 #Myr

    # Is it data?
    abstract type Data{T} end
    # Generic methods to allow broadcasting and comparison
    Base.length(x::Data) = 1
    Base.iterate(x::Data) = (x, nothing)
    Base.iterate(x::Data, state) = nothing
    Base.:(==)(x::Data, y::Data) = false
    function Base.:(==)(x::T, y::T) where {T<:Data}
        for n in fieldnames(T)
            isequal(getfield(x, n), getfield(y, n)) || return false
        end
        return true
    end

    # Reduced data
    include("analysis.jl")
    export value, stdev
    export Analysis, Analysis1D, Analysis2D, Analysis3D, Analysis4D, Analysis5D
    export CI, Age, Interval, Ellipse

    # Fitting and interpreting data
    include("regression.jl")
    export chauvenet
    export wmean, awmean, gwmean, distwmean, mswd
    export lsqfit, yorkfit, YorkFit

    include("isochron.jl")
    export age, ratio, isochron

    include("U-Pb.jl")
    export UPbAnalysis, discordance, age68, age75, stacey_kramers

    include("concordia.jl")
    export upperintercept, lowerintercept, intercepts, isconcordant, age76, ageconcordia

    include("U-Th.jl")
    include("Re-Os.jl")
    include("Lu-Hf.jl")
    include("Sm-Nd.jl")
    include("Rb-Sr.jl")
    include("K-Ar.jl")
    export UThAnalysis, ReOsAnalysis, LuHfAnalysis, SmNdAnalysis, RbSrAnalysis

    # Raw data
    include("datareduction.jl")
    export importsimsdata, calibration, calibrate, calibrate_blockwise

    include("generic_plotting.jl")
    export concordiacurve, concordiacurve!, concordialine, concordialine!,
    rankorder, rankorder!

    include("metropolis.jl")
    export metropolis_min, metropolis_min!,
    metropolis_minmax, metropolis_minmax!

    include("distributions.jl")
    export UniformDistribution, TriangularDistribution,
    HalfNormalDistribution, ExponentialDistribution,
    MeltsZirconDistribution, MeltsVolcanicZirconDistribution

    include("show.jl")

    # extra exports for pkg extensions
    export Collection, Data, AbstractAnalysis, datalimits

end # module Isoplot
