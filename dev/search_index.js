var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Isoplot","category":"page"},{"location":"#Isoplot","page":"Home","title":"Isoplot","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Isoplot.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Isoplot]","category":"page"},{"location":"#Isoplot.UPbAnalysis","page":"Home","title":"Isoplot.UPbAnalysis","text":"struct UPbAnalysis{T} <: Analysis{T}\n\nCore type for U-Pb analyses. Has fields\n\nμ :: Vector{T<:AbstractFloat}\nσ :: Vector{T<:AbstractFloat}\nΣ :: Matrix{T<:AbstractFloat}\n\nwhere μ contains the means\n\nμ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]\n\nwhere σ contains the standard deviations\n\nσ = [σ²⁰⁷Pb²³⁵U, σ²⁰⁶Pb²³⁸U]\n\nand Σ contains the covariance matrix\n\nΣ = [σ₇_₅^2 σ₇_₅*σ₃_₈\n     σ₇_₅*σ₃_₈ σ₃_₈^2]\n\nIf σ is not provided, it will be automatically calculated from Σ, given that σ.^2 = diag(Σ).\n\n\n\n\n\n","category":"type"},{"location":"#Isoplot.UPbAnalysis-NTuple{5, Number}","page":"Home","title":"Isoplot.UPbAnalysis","text":"UPbAnalysis(r²⁰⁷Pb²³⁵U, σ²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U, σ²⁰⁶Pb²³⁸U, correlation; T=Float64)\n\nConstruct a UPbAnalysis object from individual isotope ratios and (1-sigma!) uncertainties.\n\nExamples\n\njulia> UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183)\nUPbAnalysis{Float64}([22.6602, 0.40864], [0.00030625000000000004 2.4746942500000003e-6; 2.4746942500000003e-6 2.8900000000000004e-8])\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.dist_ll-Tuple{Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, Number, Number}","page":"Home","title":"Isoplot.dist_ll","text":"dist_ll(dist::Collection, mu::Collection, sigma::Collection, tmin::Number, tmax::Number)\ndist_ll(dist::Collection, analyses::Collection{<:Measurement}, tmin::Number, tmax::Number)\n\nReturn the log-likelihood of a set of mineral ages with means mu and uncertianty sigma being drawn from a given source (i.e., crystallization / closure) distribution dist, with terms to prevent runaway at low N.\n\nExamples\n\nmu, sigma = collect(100:0.1:101), 0.01*ones(11)\nll = dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100, 101)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.lsqfit-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{<:Number}}} where T<:Number","page":"Home","title":"Isoplot.lsqfit","text":"(a,b) = lsqfit(x::AbstractVector, y::AbstractVector)\n\nReturns the coefficients for a simple linear least-squares regression of the form y = a + bx\n\nExamples\n\njulia> a, b = lsqfit(1:10, 1:10)\n2-element Vector{Float64}:\n -1.19542133983862e-15\n  1.0\n\njulia> isapprox(a, 0, atol = 1e-12)\ntrue\n\njulia> isapprox(b, 1, atol = 1e-12)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.metropolis_min!-Tuple{DenseArray, Integer, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, AbstractArray, AbstractArray}","page":"Home","title":"Isoplot.metropolis_min!","text":"metropolis_min!(tminDist::DenseArray, nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)\nmetropolis_min!(tmindist::DenseArray, t0dist::DenseArray, nsteps::Integer, dist::Collection, analyses::Collection{<:UPbAnalysis}; burnin::Integer = 0) where {T}\n\nIn-place (non-allocating) version of metropolis_min, fills existing array tminDist.\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\nmetropolis_min!(tminDist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.metropolis_min-Tuple{Integer, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, Union{Tuple{Vararg{var\"#s26\", N}}, AbstractArray{var\"#s26\"}} where {var\"#s26\"<:Measurements.Measurement, N}}","page":"Home","title":"Isoplot.metropolis_min","text":"metropolis_min(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; burnin::Integer=0)\nmetropolis_min(nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)\nmetropolis_min(nsteps::Integer, dist::Collection, analyses::Collection{<:UPbAnalysis; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the minimum of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\ntmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\ntmindist, t0dist = metropolis_min(2*10^5, HalfNormalDistribution, analyses, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.metropolis_minmax!-Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, Integer, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, AbstractArray, AbstractArray}","page":"Home","title":"Isoplot.metropolis_minmax!","text":"metropolis_minmax!(tminDist, tmaxDist, llDist, acceptanceDist, nsteps::Integer, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nIn-place (non-allocating) version of metropolis_minmax, filling existing arrays\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\nmetropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.metropolis_minmax-Tuple{Integer, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}, Union{Tuple{Vararg{var\"#s26\", N}}, AbstractArray{var\"#s26\"}} where {var\"#s26\"<:Measurements.Measurement, N}}","page":"Home","title":"Isoplot.metropolis_minmax","text":"metropolis_minmax(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; burnin::Integer=0)\nmetropolis_minmax(nsteps::Integer, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)\n\nRun a Metropolis sampler to estimate the extrema of a finite-range source distribution dist using samples drawn from that distribution – e.g., estimate zircon saturation and eruption ages from a distribution of zircon crystallization ages.\n\nExamples\n\ntmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.mswd-Union{Tuple{T}, Tuple{Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where N, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}}} where T","page":"Home","title":"Isoplot.mswd","text":"mswd(μ, σ)\nmswd(μ ± σ)\n\nReturn the Mean Square of Weighted Deviates (AKA the reduced chi-squared statistic) of a dataset with values x and one-sigma uncertainties σ\n\nExamples\n\njulia> x = randn(10)\n10-element Vector{Float64}:\n -0.977227094347237\n  2.605603343967434\n -0.6869683962845955\n -1.0435377148872693\n -1.0171093080088411\n  0.12776158554629713\n -0.7298235147864734\n -0.3164914095249262\n -1.44052961622873\n  0.5515207382660242\n\njulia> mswd(x, ones(10))\n1.3901517474017941\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.wmean-Union{Tuple{T}, Tuple{Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where N, Union{Tuple{Vararg{T, N}}, AbstractArray{T}} where {T, N}}} where T","page":"Home","title":"Isoplot.wmean","text":"wμ, wσ, mswd = wmean(μ, σ; corrected=false)\nwμ ± wσ, mswd = wmean(μ ± σ; corrected=false)\n\nThe weighted mean, with or without the \"geochronologist's MSWD correction\" to uncertainty. You may specify your means and standard deviations either as separate vectors μ and σ, or as a single vector x of Measurements equivalent to x = μ .± σ\n\nIn all cases, σ is assumed to reported as actual sigma (i.e., 1-sigma).\n\nIf corrected=true, the resulting uncertainty of the weighted mean is corrected for dispersion when the MSWD is greater than 1 by multiplying by the square root of the MSWD.\n\nExamples\n\njulia> x = randn(10)\n10-element Vector{Float64}:\n  0.4612989881720301\n -0.7255529837975242\n -0.18473979056481055\n -0.4176427262202118\n -0.21975911391551833\n -1.6250003193791873\n -1.6185557291787287\n  0.25315988825847513\n -0.4979804844182867\n  1.3565281078086726\n\njulia> y = ones(10);\n\njulia> wmean(x, y)\n(-0.321824416323509, 0.31622776601683794, 0.8192171477885678)\n\njulia> wmean(x .± y)\n(-0.32 ± 0.32, 0.8192171477885678)\n\njulia> wmean(x .± y./10)\n(-0.322 ± 0.032, 81.9217147788568)\n\njulia> wmean(x .± y./10, corrected=true)\n(-0.32 ± 0.29, 81.9217147788568)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.yorkfit-Tuple{Vector{<:Measurements.Measurement}, Vector{<:Measurements.Measurement}}","page":"Home","title":"Isoplot.yorkfit","text":"yorkfit(x, σx, y, σy)\n\nUses the York (1968) two-dimensional least-squares fit to calculate a, b, and uncertanties σa, σb for the equation y = a + bx, given x, y and uncertaintes σx, `σy.\n\nFor further reference, see: York, Derek (1968) \"Least squares fitting of a straight line with correlated errors\" Earth and Planetary Science Letters 5, 320-324. doi: 10.1016/S0012-821X(68)80059-7\n\nExamples\n\njulia> x = (1:100) .+ randn.();\n\njulia> y = 2*(1:100) .+ randn.();\n\njulia> yorkfit(x, ones(100), y, ones(100))\nYorkFit{Float64}:\nLeast-squares linear fit of the form y = a + bx where\n  intercept a : -0.29 ± 0.2 (1σ)\n  slope b     : 2.0072 ± 0.0035 (1σ)\n  MSWD        : 0.8136665223891004\n\n\n\n\n\n","category":"method"}]
}
