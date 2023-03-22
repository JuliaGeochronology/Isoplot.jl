module Isoplot

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
    export ellipse

    include("regression.jl")
    export wmean, awmean, gwmean, mswd
    export lsqfit, yorkfit

    include("upb.jl")
    export UPbAnalysis, age, discordance

    include("concordia.jl")
    export upperintercept, lowerintercept, intercepts

    include("plotting.jl")
    export concordiacurve!

    include("metropolis.jl")
    export metropolis_min, metropolis_min!

    include("show.jl")

end # module Isoplot
