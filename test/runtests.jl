module BaseTests

    using Test, Statistics
    using Measurements
    using Isoplot

    @testset "Show" begin
        yf = Isoplot.YorkFit(1±1, 1±1, 0.0, 1±1, 1.0)
        @test display(yf) != NaN

        ci = CI(1:10)
        @test ci == CI{Float64}(5.5, 3.0276503540974917, 5.5, 1.225, 9.775)
        @test "$ci" === "5.5 +4.3/-4.3"
        @test display(ci) != NaN

        ci = CI(randn(100))
        @test display(ci) != NaN
    end

    @testset "General" begin
        @test Age(0) isa Age{Float64}
        @test Age(0, 1) isa Age{Float64}
        @test Interval(0, 1) isa Interval{Float64}
        @test Interval(0, 1) === Interval(0,0, 1,0)
        @test min(Interval(0, 1)) isa Age{Float64}
        @test max(Interval(0, 1)) isa Age{Float64}
        @test min(Interval(0, 1)) === Age(0,0) === Age(0)
        @test max(Interval(0, 1)) === Age(1,0) === Age(1)
        @test Isoplot.val(1) === 1
        @test Isoplot.err(1) === 0
        @test Isoplot.val(1±1) === 1.0
        @test Isoplot.err(1±1) === 1.0
        ci = CI(1:10)
        @test Isoplot.val(ci) ≈ 5.5
        @test Isoplot.err(ci) ≈ 3.0276503540974917
        a = Age(ci)
        @test Isoplot.val(a) ≈ 5.5
        @test Isoplot.err(a) ≈ 3.0276503540974917
    end

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
        @test !isnan(d1)

        a75, a68 = age(d1)
        @test a75.val ≈ 3209.725483265418
        @test a75.err ≈ 1.9420875256761048
        @test a68.val ≈ 2208.7076248184017
        @test a68.err ≈ 1.422824131349332
        @test age68(d1) == a68
        @test age75(d1) == a75
        @test discordance(d1) ≈ 31.187024051310136

        @test rand(d1) isa Vector{Float64}
        @test rand(d1, 10) isa Matrix{Float64}
        @test rand(d1, 5, 5) isa Matrix{Vector{Float64}}

        x = [22.70307499779583, 22.681635852743682, 22.63876085494785, 22.61732500220417, 22.638764147256317, 22.68163914505215, 22.70307499779583]
        y = [0.4089925091091875, 0.40901969166358015, 0.4086701825543926, 0.40829349089081246, 0.4082663083364198, 0.4086158174456074, 0.4089925091091875]
        e1 = Ellipse(d1, npoints=7)
        @test e1 isa Isoplot.Ellipse
        @test e1.x ≈ x
        @test e1.y ≈ y

        tₗₗ = 35
        ui = upperintercept(tₗₗ ± 10, d1)
        @test ui == 3921.343026090256 ± 2.745595368456398
        ui = upperintercept(tₗₗ, d1)
        @test ui == 3921.343026090256 ± 0.7241111646504936

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


        # Stacey-Kramers common Pb model
        @test stacey_kramers(0) == (18.7, 15.628)
        @test stacey_kramers(3700) == (11.152, 12.998)
        @test stacey_kramers(4567) == (9.314476625036953, 12.984667029161916)
        @test stacey_kramers(5000) === (NaN, NaN)

    end
    @testset "Other systems" begin
        μ, σ = rand(2), rand(2)
        @test UThAnalysis(μ, σ) isa UThAnalysis
        @test ReOsAnalysis(μ, σ) isa ReOsAnalysis
        @test LuHfAnalysis(μ, σ) isa LuHfAnalysis
        @test SmNdAnalysis(μ, σ) isa SmNdAnalysis
        @test RbSrAnalysis(μ, σ) isa RbSrAnalysis
    end

    @testset "Weighted means" begin
        # Weighted means
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

        N = 10^6
        c = randn(N).+50
        a,b = randn(N), randn(N).+1
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,1], [1,1]; corrected=false)
        @test mean(d) ≈ μ atol = 0.02
        @test std(d) ≈ σ atol = 0.002
        d = distwmean(a,b; corrected=true)
        @test mean(d) ≈ μ atol = 0.1
        @test std(d) ≈ std(vcat(a,b)) atol = 0.005

        b .+= 9
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,10], [1,1]; corrected=false)
        @test mean(d) ≈ μ atol = 0.02
        @test std(d) ≈ σ atol = 0.002
        d = distwmean(a,b; corrected=true)
        μ,σ,_ = wmean([0,10], [1,1]; corrected=true)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.2
        @test std(d) ≈ std(vcat(a,b)) atol = 0.02

        b .+= 90
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,100], [1,1]; corrected=false)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.02
        d = distwmean(a,b; corrected=true)
        μ,σ,_ = wmean([0,100], [1,1]; corrected=true)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.2
        @test std(d) ≈ std(vcat(a,b)) atol = 0.2

        d = distwmean(a,b,c; corrected=false)
        μ,σ,_ = wmean([0,50,100], [1,1,1]; corrected=false)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.02
        d = distwmean(a,b,c; corrected=true)
        μ,σ,_ = wmean([0,50,100], [1,1,1]; corrected=true)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 2


        a,b = 2randn(N), 3randn(N).+1
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,1], [2,3]; corrected=false)
        @test mean(d) ≈ μ atol = 0.05
        @test std(d) ≈ σ atol = 0.005

        b .+= 9
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,10], [2,3]; corrected=false)
        @test mean(d) ≈ μ atol = 0.05
        @test std(d) ≈ σ atol = 0.005
        d = distwmean(a,b; corrected=true)
        μ,σ,_ = wmean([0,10], [2,3]; corrected=true)
        @test mean(d) ≈ μ atol = 0.1
        @test std(d) ≈ σ atol = 1

        b .+= 90
        d = distwmean(a,b; corrected=false)
        μ,σ,_ = wmean([0,100], [2,3]; corrected=false)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.02
        d = distwmean(a,b; corrected=true)
        μ,σ,_ = wmean([0,100], [2,3]; corrected=true)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 2

        d = distwmean(a,b,c; corrected=false)
        μ,σ,_ = wmean([0,50,100], [2,1,3]; corrected=false)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.02
        d = distwmean(a,b,c; corrected=true)
        μ,σ,_ = wmean([0,50,100], [2,1,3]; corrected=true)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 2
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

    @testset "Regression" begin
        # Simple linear regression
        ϕ = lsqfit(1:10, 1:10)
        @test ϕ[1] ≈ 0 atol=1e-12
        @test ϕ[2] ≈ 1 atol=1e-12

        # York (1968) fit
        x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
        y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
        σx = σy = ones(10)/4
        yf = yorkfit(x, σx, y, σy)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈-0.23498964673701916
        @test yf.intercept.err ≈ 0.02250863813481163
        @test yf.slope.val ≈ 1.041124018512526
        @test yf.slope.err ≈ 0.0035683808205783673
        @test yf.mswd ≈ 1.1419901440278089

        yf = yorkfit(x, σx, y, σy, zeros(length(x)))
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈-0.2446693790977319
        @test yf.intercept.err ≈ 0.2469541320914601
        @test yf.slope.val ≈ 1.0428730084538775
        @test yf.slope.err ≈ 0.039561084436542084
        @test yf.mswd ≈ 1.1417951538670306

        x = ((1:100) .+ randn.()) .± 1
        y = (2*(1:100) .+ randn.()) .± 1
        yf = yorkfit(x, y)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈ 0 atol = 2
        @test yf.slope.val ≈ 2 atol = 0.1
        @test yf.mswd ≈ 1 atol = 0.5

        yf = yorkfit(analyses)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈ 0.0050701916562521515
        @test yf.intercept.err ≈ 0.004099648408656529
        @test yf.slope.val ≈ 0.1079872513087868
        @test yf.slope.err ≈ 0.0037392146940848233
        @test yf.xm ≈ 1.096376584184683
        @test yf.ym.val ≈ 0.12346488538167279
        @test yf.ym.err ≈ 2.235949353726133e-5
        @test yf.mswd ≈ 0.41413597765872123

        ui = upperintercept(analyses)
        @test ui.val ≈ 752.6744316220871
        @test ui.err ≈ 0.5288009504134864

        li = lowerintercept(analyses)
        @test li.val ≈ 115.83450556482211
        @test li.err ≈ 94.4384248140631

        ui, li = intercepts(analyses)
        @test ui.val ≈ 752.6744316220871
        @test ui.err ≈ 0.5288009504134864
        @test li.val ≈ 115.83450556482211
        @test li.err ≈ 94.4384248140631

    end

    @testset "Concordia Metropolis" begin
        data = upperintercept.(0, analyses)
        @test Isoplot.dist_ll(ones(10), data, 751, 755) ≈ -20.09136536048026
        @test Isoplot.dist_ll(ones(10), data, 750, 760) ≈ -30.459633175497830
        @test Isoplot.dist_ll(ones(10), data, 752, 753) ≈ -15.305463167234748
        @test Isoplot.dist_ll(ones(10), data, 751, 752) ≈ -47.386667785224034

        tmindist, t0dist = metropolis_min(1000, ones(10), analyses; burnin=200)
        @test tmindist isa Vector{Float64}
        @test mean(tmindist) ≈ 751.85 atol = 1.5
        @test std(tmindist) ≈ 0.40 rtol = 0.6
        @test t0dist isa Vector{Float64}
        @test mean(t0dist) ≈ 80. atol = 90
        @test std(t0dist) ≈ 50. rtol = 0.6

        tmindist, tmaxdist, t0dist, lldist, acceptancedist = metropolis_minmax(10000, ones(10), analyses; burnin=200)
        @test tmindist isa Vector{Float64}
        @test mean(tmindist) ≈ 751.85 atol = 1.5
        @test std(tmindist) ≈ 0.40 rtol = 0.6
        @test tmaxdist isa Vector{Float64}
        @test mean(tmaxdist) ≈ 753.32 atol = 1.5
        @test std(tmaxdist) ≈ 0.60 rtol = 0.6
        @test t0dist isa Vector{Float64}
        @test mean(t0dist) ≈ 80. atol = 90
        @test std(t0dist) ≈ 50. rtol = 0.6
        @test lldist isa Vector{Float64}
        @test acceptancedist isa BitVector
        @test mean(acceptancedist) ≈ 0.6 atol=0.2

        terupt = CI(tmindist)
        @test terupt isa CI{Float64}
        @test terupt.mean ≈ 751.85 atol = 1.5
        @test terupt.sigma ≈ 0.40 rtol = 0.6
        @test terupt.median ≈ 751.83 atol = 1.5
        @test terupt.lower ≈ 750.56 atol = 1.5
        @test terupt.upper ≈ 752.52 atol = 1.5
    end

    @testset "General Metropolis" begin
        mu, sigma = collect(100:0.1:101), 0.01*ones(11);
        @test Isoplot.dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100,101) ≈ -3.6933372932657607

        tmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu .± sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228 atol=0.015

        tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu .± sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test mean(tmaxdist) ≈ 101.08  atol=0.015
        @test lldist isa Vector{Float64}
        @test acceptancedist isa BitVector
        @test mean(acceptancedist) ≈ 0.6 atol=0.2

        @test mean(UniformDistribution) ≈ 1
        @test mean(TriangularDistribution) ≈ 1
        @test mean(HalfNormalDistribution) ≈ 1
        @test mean(ExponentialDistribution) ≈ 1.03 atol=0.01
        @test mean(MeltsZirconDistribution) ≈ 1 atol=0.01
        @test mean(MeltsVolcanicZirconDistribution) ≈ 1 atol=0.01
    end

end

module PlotsTest

    using Test, Statistics
    using Measurements
    using Isoplot
    using ImageIO, FileIO,Plots

    import ..BaseTests: analyses

    # Base.retry_load_extensions()
    @testset "Plotting" begin
        # Plot single concordia ellipse
        h = plot(analyses[1], color=:blue, alpha=0.3, label="", framestyle=:box)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.8151617156862782,0.8151617156862782,0.986395212418301) rtol = 0.02
        rm("concordia.png")

        # Plot many concordia ellipses and concordia curve
        h = plot(analyses, color=:blue, alpha=0.3, label="", framestyle=:box)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9448600653594766,0.9448600653594766,0.9658495915032675) rtol = 0.01
        rm("concordia.png")

        # Plot single concordia line
        h = concordialine(0, 100, label="")
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.981812385620915,0.9832414705882352,0.9841176307189541) rtol = 0.01
        rm("concordia.png")
        h = concordialine(0, 100, label="", truncate=true)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.981812385620915,0.9832414705882352,0.9841176307189541) rtol = 0.01
        rm("concordia.png")

        # Plot many single concordia lines
        h = concordialine(10*randn(100).+10, 100*randn(100).+1000, label="")
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.966250588235294,0.966250588235294,0.966250588235294) rtol = 0.01
        rm("concordia.png")
        h = concordialine(10*randn(100).+10, 100*randn(100).+1000, label="", truncate=true)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.966250588235294,0.966250588235294,0.966250588235294) rtol = 0.01
        rm("concordia.png")

        # Rank-order plot
        h = rankorder(1:10, 2*ones(10), label="")
        savefig(h, "rankorder.png")
        img = load("rankorder.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9842757843137254,0.9882103758169932,0.9906218464052285) rtol = 0.01
        rm("rankorder.png")
    end

end

module MakieTest
    using Test, Statistics
    using Measurements
    using Isoplot
    using ImageIO, FileIO, CairoMakie
    using ColorTypes

    import ..BaseTests: analyses

    # Base.retry_load_extensions()
    @testset "Makie Plotting" begin
        f = Figure()
        ax = Axis(f[1,1])
        plot!(analyses[1], color=(:blue,0.3))
        save("concordia.png",f)
        img = load("concordia.png")
        @test size(img) == (900, 1200)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.8524913580246877,0.8524913580246877,0.9885884168482209) rtol = 0.02
        rm("concordia.png")
    
        # Plot many concordia ellipses and concordia curve
        f2 = Figure()
        ax2 = Axis(f2[1,1])
        plot!.(analyses, color=(:blue, 0.3))
        ages = age.(analyses)
        concordiacurve!(minimum(ages)[1].val-5,maximum(ages)[1].val+5)
        
        xmin, xmax, ymin, ymax = datalimits(analyses)
        limits!(ax2,xmin,xmax,ymin,ymax)
        save("concordia.png",f2)
        img = load("concordia.png")
        @test size(img) == (900, 1200)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9523360547065816,0.9523360547065816,0.9661779080315414) rtol = 0.01
        rm("concordia.png")
    
        # Plot single concordia line
        f3 = Figure()
        ax3 = Axis(f3[1,1])
        concordialine!(0, 100)
        concordiacurve!(0,100)
        save("concordia.png",f3)
        img = load("concordia.png")
        @test size(img) == (900, 1200)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9845678970995279,0.9845678970995279,0.9845678970995279) rtol = 0.01
        rm("concordia.png")
       
    end
end
