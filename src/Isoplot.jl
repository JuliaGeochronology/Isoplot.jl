module Isoplot

    using VectorizedStatistics
    using LoopVectorization
    using LinearAlgebra
    using Distributions
    using Measurements
    using Plots: Shape, plot, plot!
    import Plots
    export plot, plot!

    # A type alias for array-ish types
    const Collection{T} = Union{AbstractArray{T}, NTuple{N,T}} where N

    # Age of Earth and the Solar System
    const t🜨 = 4.567e3 #Myr

    # Abstract types which we'll subtype later
    include("analysis.jl")
    export age, ratio, ellipse, CI, val, err

    include("regression.jl")
    export wmean, awmean, gwmean, distwmean, mswd
    export lsqfit, yorkfit

    include("U-Pb.jl")
    export UPbAnalysis, discordance, age68, age75

    include("concordia.jl")
    export upperintercept, lowerintercept, intercepts

    include("U-Th.jl")
    include("Re-Os.jl")
    include("Lu-Hf.jl")
    include("Sm-Nd.jl")
    include("Rb-Sr.jl")
    include("K-Ar.jl")
    export UThAnalysis, ReOsAnalysis, LuHfAnalysis, SmNdAnalysis, RbSrAnalysis

    include("plotting.jl")
    export concordiacurve!, concordialine, concordialine!,
    rankorder, rankorder!

    include("metropolis.jl")
    export metropolis_min, metropolis_min!,
    metropolis_minmax, metropolis_minmax!

    include("distributions.jl")
    export UniformDistribution, TriangularDistribution,
    HalfNormalDistribution, ExponentialDistribution,
    MeltsZirconDistribution, MeltsVolcanicZirconDistribution

    include("show.jl")

end # module Isoplot
