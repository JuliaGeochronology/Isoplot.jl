var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Isoplot","category":"page"},{"location":"#Isoplot","page":"Home","title":"Isoplot","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Isoplot.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Isoplot]","category":"page"},{"location":"#Isoplot.UPbAnalysis","page":"Home","title":"Isoplot.UPbAnalysis","text":"struct UPbAnalysis{T} <: Analysis{T}\n\nCore type for U-Pb analyses. Has fields\n\nμ :: Vector{T<:AbstractFloat}\nΣ :: Matrix{T<:AbstractFloat}\n\nwhere μ contains the means\n\nμ = [r²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U]\n\nand Σ contains the covariance matrix\n\nΣ = [σ₇_₅^2 σ₇_₅*σ₃_₈\n     σ₇_₅*σ₃_₈ σ₃_₈^2]\n\n\n\n\n\n","category":"type"},{"location":"#Isoplot.UPbAnalysis-NTuple{5, Number}","page":"Home","title":"Isoplot.UPbAnalysis","text":"UPbAnalysis(r²⁰⁷Pb²³⁵U, σ²⁰⁷Pb²³⁵U, r²⁰⁶Pb²³⁸U, σ²⁰⁶Pb²³⁸U, correlation; T=Float64)\n\nConstruct a UPbAnalysis object from individual isotope ratios and (1-sigma!) uncertainties.\n\nExamples\n\njulia> UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183)\nUPbAnalysis{Float64}([22.6602, 0.40864], [0.00030625000000000004 2.4746942500000003e-6; 2.4746942500000003e-6 2.8900000000000004e-8])\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.wmean-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractArray}} where T","page":"Home","title":"Isoplot.wmean","text":"wμ, wσ, mswd = wmean(μ, σ; corrected=false)\nwμ ± wσ, mswd = wmean(x; corrected=false)\n\nThe weighted mean, with or without the \"geochronologist's MSWD correction\" to uncertainty. You may specify your means and standard deviations either as separate vectors μ and σ, or as a single vector x of Measurements equivalent to x = μ .± σ\n\nIn all cases, σ is assumed to reported as actual sigma (i.e., 1-sigma).\n\nExamples\n\njulia> x = randn(10)\n10-element Vector{Float64}:\n  0.4612989881720301\n -0.7255529837975242\n -0.18473979056481055\n -0.4176427262202118\n -0.21975911391551833\n -1.6250003193791873\n -1.6185557291787287\n  0.25315988825847513\n -0.4979804844182867\n  1.3565281078086726\n\njulia> y = ones(10);\n\njulia> wmean(x, y)\n(-0.321824416323509, 0.31622776601683794, 0.8192171477885678)\n\njulia> wmean(x .± y)\n(-0.32 ± 0.32, 0.8192171477885678)\n\njulia> wmean(x .± y./10)\n(-0.322 ± 0.032, 81.9217147788568)\n\njulia> wmean(x .± y./10, corrected=true)\n(-0.32 ± 0.29, 81.9217147788568)\n\n\n\n\n\n","category":"method"},{"location":"#Isoplot.yorkfit-Tuple{Vector{<:Measurements.Measurement}, Vector{<:Measurements.Measurement}}","page":"Home","title":"Isoplot.yorkfit","text":"yorkfit(x, σx, y, σy)\n\nUses the York (1968) two-dimensional least-squares fit to calculate a, b, and uncertanties σa, σb for the equation y = a + bx, given x, y and uncertaintes σx, `σy.\n\nFor further reference, see: York, Derek (1968) \"Least squares fitting of a straight line with correlated errors\" Earth and Planetary Science Letters 5, 320-324. doi: 10.1016/S0012-821X(68)80059-7\n\nExamples\n\njulia> x = (1:100) .+ randn.();\n\njulia> y = 2*(1:100) .+ randn.();\n\njulia> yorkfit(x, ones(100), y, ones(100))\nYorkFit{Float64}:\nLeast-squares linear fit of the form y = a + bx where\n  intercept a : -0.29 ± 0.2 (1σ)\n  slope b     : 2.0072 ± 0.0035 (1σ)\n  MSWD        : 0.8136665223891004\n\n\n\n\n\n","category":"method"}]
}
