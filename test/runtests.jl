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
    @test d1 isa UPbAnalysis{Float64}
    @test d2 isa UPbAnalysis{Float64}
    @test d1.μ ≈ d2.μ
    @test d1.Σ ≈ d2.Σ

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

    # Simple linear egression
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
    @test fobj.intercept.val ≈ 0 atol = 1
    @test fobj.slope.val ≈ 2 atol = 0.1
    @test fobj.mswd ≈ 1 atol = 0.5

end

using ImageIO, FileIO
@testset "Plotting" begin
    data = [ 1.02849 0.0006788034 0.11781 3.946635e-5 0.959
             1.02826 0.0013932923 0.11771 5.120385e-5 0.693
             1.02716 0.0006727898 0.11768 3.8246e-5 0.972
             1.02694 0.001899839 0.11767 9.707775e-5 0.65
             1.02712 0.0009090012 0.11763 4.881645e-5 0.797
             1.02551 0.00081528045 0.11756 4.1146e-5 0.859
             1.02653 0.0014679379 0.11754 7.69887e-5 0.691
             1.02451 0.0006556864 0.11742 3.75744e-5 0.962
             1.02039 0.00309688365 0.11738 9.21433e-5 0.552
             1.01931 0.0008460273 0.11699 4.21164e-5 0.837
             1.01996 0.000662974 0.11689 3.798925e-5 0.971
             1.01914 0.0007235894 0.11689 3.798925e-5 0.929
             1.01309 0.0006686394 0.11614 3.77455e-5 0.962
             1.01285 0.001397733 0.11594 5.33324e-5 0.645
             1.00183 0.0008315189 0.11502 3.91068e-5 0.856
             1.00198 0.0006663167 0.11496 3.7362e-5 0.968  ]

    d = UPbAnalysis.(eachcol(data)...,)

    # Plot many concordia ellipses and concordia curve
    h = plot(framestyle=:box)
    plot!(h, d, color=:blue, alpha=0.3, label="")
    concordiacurve!(h)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9448600653594766,0.9448600653594766,0.9658495915032675) rtol = 0.01
    rm("concordia.png")

    h = plot(d, color=:blue, alpha=0.3, label="", framestyle=:box)
    concordiacurve!(h)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9448600653594766,0.9448600653594766,0.9658495915032675) rtol = 0.01
    rm("concordia.png")

    # Plot single concordia ellipse
    h = plot(framestyle=:box)
    plot!(h, d[1], color=:blue, alpha=0.3, label="")
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9414076797385613,0.9414076797385613,0.9870670098039216) rtol = 0.01
    rm("concordia.png")

    h = plot(d[1], color=:blue, alpha=0.3, label="", framestyle=:box)
    savefig(h, "concordia.png")
    img = load("concordia.png")
    @test size(img) == (400,600)
    @test sum(img)/length(img) ≈ RGB{Float64}(0.9414076797385613,0.9414076797385613,0.9870670098039216) rtol = 0.01
    rm("concordia.png")
end
