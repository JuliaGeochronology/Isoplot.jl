module Isoplot

    using VectorizedStatistics
    using LoopVectorization: @turbo
    using LinearAlgebra
    using Distributions
    using Measurements
    

    # A type alias for array-ish types
    const Collection{T} = Union{AbstractArray{T}, NTuple{N,T}} where N

    # Age of Earth and the Solar System
    const t🜨 = 4.567e3 #Myr

    # Abstract types which we'll subtype later
    include("analysis.jl")
    export age, ratio, CI, Age, Interval, Ellipse

    include("regression.jl")
    export wmean, awmean, gwmean, distwmean, mswd
    export lsqfit, yorkfit

    include("U-Pb.jl")
    export UPbAnalysis, discordance, age68, age75, stacey_kramers

    include("concordia.jl")
    export upperintercept, lowerintercept, intercepts

    include("U-Th.jl")
    include("Re-Os.jl")
    include("Lu-Hf.jl")
    include("Sm-Nd.jl")
    include("Rb-Sr.jl")
    include("K-Ar.jl")
    export UThAnalysis, ReOsAnalysis, LuHfAnalysis, SmNdAnalysis, RbSrAnalysis

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

    #extra exports for pkg extenstions
    export Data, Analysis, Collection, val, err, vminimum, vmaximum, datalimits

end # module Isoplot
