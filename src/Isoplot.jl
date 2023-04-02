module Isoplot

    using VectorizedStatistics
    using LoopVectorization
    using LinearAlgebra
    using Distributions
    using Measurements
    using Plots: Shape, plot, plot!
    import Plots
    export plot, plot!

    const Collection{T} = Union{AbstractArray{T}, NTuple{N,T}} where N

    # Abstract types which we'll subtype later
    include("analysis.jl")
    export age, ratio, ellipse, CI

    include("regression.jl")
    export wmean, awmean, gwmean, mswd
    export lsqfit, yorkfit

    include("U-Pb.jl")
    export UPbAnalysis, discordance

    include("concordia.jl")
    export upperintercept, lowerintercept, intercepts

    include("U-Th.jl")
    include("Re-Os.jl")
    include("Lu-Hf.jl")
    include("Sm-Nd.jl")
    include("Rb-Sr.jl")
    include("K-Ar.jl")

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
