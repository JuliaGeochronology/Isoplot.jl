using Isoplot
using Test, Statistics
using Plots
using Measurements

@testset "U-Pb" begin
    r75 = 22.6602
    σ75 = 0.017516107998
    r68 = 0.408643
    σ68 = 0.0001716486532565
    corr = 0.831838

    d1 = UPbAnalysis(r75, σ75, r68, σ68, corr)
    d2 = UPbAnalysis([22.6602, 0.408643], [0.00030681403939759964 2.501017729814154e-6; 2.501017729814154e-6 2.9463260164770177e-8])
    d3 = UPbAnalysis([22.6602, 0.408643], [0.017516107998, 0.0001716486532565], [0.00030681403939759964 2.501017729814154e-6; 2.501017729814154e-6 2.9463260164770177e-8])
    @test d1 isa UPbAnalysis{Float64}
    @test d2 isa UPbAnalysis{Float64}
    @test d3 isa UPbAnalysis{Float64}
    @test d1.μ ≈ d2.μ ≈ d3.μ
    @test d1.σ ≈ d2.σ ≈ d3.σ
    @test d1.Σ ≈ d2.Σ ≈ d3.Σ

    a76, a68 = age(d1)
    @test a76.val ≈ 3209.725483265418
    @test a76.err ≈ 0.933031271855386
    @test a68.val ≈ 2208.7076248184017
    @test a68.err ≈ 1.422824131349332
    @test discordance(d1) ≈ 31.187024051310136

    @test rand(d1) isa Vector{Float64}
    @test rand(d1, 10) isa Matrix{Float64}

    x = [22.70307499779583, 22.681635852743682, 22.63876085494785, 22.61732500220417, 22.638764147256317, 22.68163914505215, 22.70307499779583]
    y = [0.4089925091091875, 0.40901969166358015, 0.4086701825543926, 0.40829349089081246, 0.4082663083364198, 0.4086158174456074, 0.4089925091091875]
    e1 = ellipse(d1, npoints=7)
    @test e1 isa Isoplot.Ellipse
    @test e1.x ≈ x
    @test e1.y ≈ y

    tₗₗ = 35
    ui = upperintercept(tₗₗ, d1)
    @test ui isa Isoplot.Measurement

    N = 10000
    uis = upperintercept(tₗₗ, d1, N)
    @test uis isa Vector{Float64}
    @test mean(uis) ≈ ui.val atol=(4*ui.err/sqrt(N))
    @test std(uis) ≈ ui.err rtol=0.03

    # Test upper and lower intercepts of multiple-sample concordia arrays
    d = [UPbAnalysis(22.6602, 0.0175, 0.40864, 0.00017, 0.83183)
         UPbAnalysis(33.6602, 0.0175, 0.50864, 0.00017, 0.83183)]

    uis = upperintercept(d, N)
    @test mean(uis) ≈ 4601.82 atol=0.1
    @test std(uis) ≈ 1.53 atol=0.1

    lis = lowerintercept(d, N)
    @test mean(lis) ≈ 1318.12 atol=0.1
    @test std(lis) ≈ 2.04 atol=0.1

    uis, lis = intercepts(d, N)
    @test mean(uis) ≈ 4601.82 atol=0.1
    @test std(uis) ≈ 1.53 atol=0.1
    @test mean(lis) ≈ 1318.12 atol=0.1
    @test std(lis) ≈ 2.04 atol=0.1
end

@testset "Regression" begin
    x = [-3.4699, -0.875, -1.4189, 1.2993, 1.1167, 0.8357, 0.9985, 1.2789, 0.5446, 0.5639]
    σx = ones(10)/4
    @test all(awmean(x, σx) .≈ (0.08737999999999996, 0.07905694150420949, 38.44179426844445))
    @test all(gwmean(x, σx) .≈ (0.08737999999999996, 0.49016447665837415, 38.44179426844445))
    @test awmean(x .± σx) == (0.08737999999999996 ± 0.07905694150420949, 38.44179426844445)
    @test gwmean(x .± σx) == (0.08737999999999996 ± 0.49016447665837415, 38.44179426844445)
    @test mswd(x, σx) ≈ 38.44179426844445
    @test mswd(x .± σx) ≈ 38.44179426844445
    σx .*= 20 # Test underdispersed data
    @test gwmean(x, σx) == awmean(x, σx)

    # Simple linear regression
    ϕ = lsqfit(1:10, 1:10)
    @test ϕ[1] ≈ 0 atol=1e-12
    @test ϕ[2] ≈ 1 atol=1e-12

    # York (1968) fit
    x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
    y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
    σx = σy = ones(10)/4
    fobj = yorkfit(x, σx, y, σy)
    @test fobj isa Isoplot.YorkFit
    @test fobj.intercept.val ≈-0.23498964673701916
    @test fobj.intercept.err ≈ 0.02250863813481163
    @test fobj.slope.val ≈ 1.041124018512526
    @test fobj.slope.err ≈ 0.0035683808205783673
    @test fobj.mswd ≈ 1.1419901440278089
    @test display(fobj) != NaN

    x = ((1:100) .+ randn.()) .± 1
    y = (2*(1:100) .+ randn.()) .± 1
    fobj = yorkfit(x, y)
    @test fobj isa Isoplot.YorkFit
    @test fobj.intercept.val ≈ 0 atol = 2
    @test fobj.slope.val ≈ 2 atol = 0.1
    @test fobj.mswd ≈ 1 atol = 0.5
end

data = [1.1009 0.00093576 0.123906 0.00002849838 0.319
        1.1003 0.00077021 0.123901 0.00003531178 0.415
        1.0995 0.00049477 0.123829 0.00002538494 0.434
        1.0992 0.00060456 0.123813 0.00003652483 0.616
        1.1006 0.00071539 0.123813 0.00002228634 0.321
        1.0998 0.00076986 0.123802 0.00002537941 0.418
        1.0992 0.00065952 0.123764 0.00003589156 0.509
        1.0981 0.00109810 0.123727 0.00003959264 0.232
        1.0973 0.00052670 0.123612 0.00002966688 0.470
        1.0985 0.00087880 0.123588 0.00002842524 0.341
        1.0936 0.00054680 0.123193 0.00003264614 0.575
        1.0814 0.00051366 0.121838 0.00003045950 0.587 ]

analyses = UPbAnalysis.(eachcol(data)...,)

using ImageIO, FileIO
@testset "Plotting" begin
    # Plot many concordia ellipses and concordia curve
    h = plot(framestyle=:box)
    plot!(h, analyses, color=:blue, alpha=0.3, label="")
    concordiacurve!(h)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9493200816993481,0.9493200816993481,0.9614482516339886) rtol = 0.01
    rm("concordia.png")

    h = plot(analyses, color=:blue, alpha=0.3, label="", framestyle=:box)
    concordiacurve!(h)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9448600653594766,0.9448600653594766,0.9658495915032675) rtol = 0.01
    rm("concordia.png")

    # Plot single concordia ellipse
    h = plot(framestyle=:box)
    plot!(h, analyses[1], color=:blue, alpha=0.3, label="")
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.8151617156862782,0.8151617156862782,0.986395212418301) rtol = 0.02
    rm("concordia.png")

    h = plot(analyses[1], color=:blue, alpha=0.3, label="", framestyle=:box)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.8151617156862782,0.8151617156862782,0.986395212418301) rtol = 0.02
    rm("concordia.png")
end

@testset "Metropolis" begin
    data = upperintercept.(0, analyses)
    @test Isoplot.dist_ll(ones(10), data, 751, 755) ≈ -20.139445921297565
    @test Isoplot.dist_ll(ones(10), data, 750, 760) ≈ -30.508892188237297
    @test Isoplot.dist_ll(ones(10), data, 752, 753) ≈ -15.347952846611049
    @test Isoplot.dist_ll(ones(10), data, 751, 752) ≈ -47.431672803923334

    tmindist, t0dist = metropolis_min(1000, ones(10), analyses; burnin=200)
    @test tmindist isa Vector{Float64}
    @test mean(tmindist) ≈ 751.8964100608977 atol = 1.5
    @test std(tmindist) ≈ 0.44203784229237997 rtol = 0.6
    @test t0dist isa Vector{Float64}
    @test mean(t0dist) ≈ 43.95580821273335 atol = 90
    @test std(t0dist) ≈ 31.489531276914047 rtol = 0.6
end
