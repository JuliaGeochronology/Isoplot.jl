module BaseTests

    using Test, Statistics
    using Measurements
    using Isoplot

    @testset "Show" begin
        yf = Isoplot.YorkFit(1±1, 1±1, 0.0, 1±1, 1.0)
        @test display(yf) != NaN

        ci = CI(1:10)
        @test ci == CI{Float64}(5.5, 3.0276503540974917, 5.5, 1.225, 9.775)
        @test "$ci" === "5.5 +4.28/-4.28"
        @test display(ci) != NaN

        ci = CI(randn(10_000))
        @test ci.mean ≈ 0 atol = 0.1
        @test ci.sigma ≈ 1 atol = 0.1
        @test ci.lower ≈ -1.9599639845400576 atol = 0.15
        @test ci.upper ≈ 1.9599639845400576 atol = 0.15
        @test display(ci) != NaN

        ci = CI([0.4227751675791465, 2.339746777712491, 0.7073822332289744, 1.115676532764763, 1.3295183295746795, 1.1839524860728812, 9.051763176763613, 3.163690972571135, 0.35529944660568286, 20.552275066037705, 1.3021911627717235, 0.8154804384732102, 2.630294225259665, 0.24662095224118538, 0.4755015782877633])
        @test ci == CI{Float64}(3.0461445697296416, 5.314497528659032,1.1839524860728812, 0.2846584252687595, 16.52709590479178)
    end

    @testset "General" begin
        # 2D analyses
        x = rand(100,2)
        d1 = Analysis(x)
        d2 = Analysis(x[:,1], x[:,2])
        d3 = Analysis(nanmean(x[:,1]), nansem(x[:,1]), nanmean(x[:,2]), nansem(x[:,2]), nancor(x[:,1],x[:,2]))
        @test d1 isa Analysis2D{Float64}
        @test d2 isa Analysis2D{Float64}
        @test d3 isa Analysis2D{Float64}
        @test d1.μ ≈ d2.μ ≈ d3.μ
        @test d1.σ ≈ d2.σ ≈ d3.σ
        @test d1.Σ ≈ d2.Σ ≈ d3.Σ
        @test !isnan(d1) && !isnan(d2) && !isnan(d3)
        # 3D analyses
        x = rand(100,3)
        d1 = Analysis(x)
        d2 = Analysis(x[:,1], x[:,2], x[:,3])
        @test d1 isa Analysis3D{Float64}
        @test d2 isa Analysis3D{Float64}
        @test d1.μ ≈ d2.μ
        @test d1.σ ≈ d2.σ
        @test d1.Σ ≈ d2.Σ
        @test !isnan(d1) && !isnan(d2)
        # Age, Interval, and CI types
        @test Age(0) isa Age{Float64}
        @test Age(0, 1) isa Age{Float64}
        @test Interval(0, 1) isa Interval{Float64}
        @test Interval(0, 1) === Interval(0,0, 1,0)
        @test min(Interval(0, 1)) isa Age{Float64}
        @test max(Interval(0, 1)) isa Age{Float64}
        @test min(Interval(0, 1)) === Age(0,0) === Age(0)
        @test max(Interval(0, 1)) === Age(1,0) === Age(1)
        @test value(1) === 1
        @test stdev(1) === 0
        @test value(1±1) === 1.0
        @test stdev(1±1) === 1.0
        ci = CI(1:10)
        @test value(ci) ≈ 5.5
        @test stdev(ci) ≈ 3.0276503540974917
        a = Age(ci)
        @test value(a) ≈ 5.5
        @test stdev(a) ≈ 3.0276503540974917
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
        @test d1.data.μ ≈ d2.data.μ ≈ d3.data.μ
        @test d1.data.σ ≈ d2.data.σ ≈ d3.data.σ
        @test d1.data.Σ ≈ d2.data.Σ ≈ d3.data.Σ
        @test !isnan(d1)

        a75, a68 = age(d1)
        @test a75.val ≈ 3209.725483265418
        @test a75.err ≈ 1.9420875256761048
        @test a68.val ≈ 2208.7076248184017
        @test a68.err ≈ 1.422824131349332
        @test age68(d1) == a68
        @test age75(d1) == a75
        @test discordance(d1) ≈ 31.187024051310136
        @test isconcordant(d1) === false

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
        @test ui.val ≈ 3921.343026090256
        @test ui.err ≈ 2.745595368456472
        ui = upperintercept(tₗₗ, d1)
        @test ui.val ≈ 3921.343026090256
        @test ui.err ≈ 0.7241111646504936

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
        x1, x2 = rand(20), rand(20)
        @test UThAnalysis(x1, x2) isa UThAnalysis
        @test ReOsAnalysis(x1, x2) isa ReOsAnalysis
        @test LuHfAnalysis(x1, x2) isa LuHfAnalysis
        @test SmNdAnalysis(x1, x2) isa SmNdAnalysis
        @test RbSrAnalysis(x1, x2) isa RbSrAnalysis
    end

    @testset "Isochrons" begin
        N = 100
        slope = 0.05

        analyses = [RbSrAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.706) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, RbSrAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=3stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.706       atol=3stdev(c.line.intercept)
        @test value(age(c)) ≈ 3502.5243481286466    atol=3stdev(age(c))
        @test stdev(age(c)) ≈ 5.104097410724326     atol=1
        @test value(age(slope, Isoplot.λ87Rb)) ≈ 3502.5243481286466

        analyses = [LuHfAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.2828) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, LuHfAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=3stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.2828      atol=3stdev(c.line.intercept)
        @test value(age(c)) ≈ 2613.2921354810956    atol=3stdev(age(c))
        @test stdev(age(c)) ≈ 11.217070579075328    atol=2
        @test value(age(slope, Isoplot.λ176Lu)) ≈ 2613.2921354810956

        analyses = [ReOsAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.14848) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, ReOsAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=3stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.14848     atol=3stdev(c.line.intercept)
        @test value(age(c)) ≈ 2923.492370389601     atol=3stdev(age(c))
        @test stdev(age(c)) ≈ 5.4768445350318045    atol=1
        @test value(age(slope, Isoplot.λ187Re)) ≈ 2923.492370389601

        analyses = [SmNdAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.5126) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, SmNdAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=3stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.5126      atol=3stdev(c.line.intercept)
        @test value(age(c)) ≈ 7478.565936454944     atol=3stdev(age(c))
        @test stdev(age(c)) ≈ 13.869217589658614    atol=3
        @test value(age(slope, Isoplot.λ147Sm)) ≈ 7478.565936454944
    end

    @testset "Weighted means" begin
        # Weighted means
        x = [-3.4699, -0.875, -1.4189, 1.2993, 1.1167, 0.8357, 0.9985, 1.2789, 0.5446, 0.5639]
        σx = ones(10)/4
        @test all(awmean(x, σx) .≈ (0.08737999999999996, 0.07905694150420949, 38.44179426844445))
        @test all(gwmean(x, σx) .≈ (0.08737999999999996, 0.49016447665837415, 38.44179426844445))
        wm, m = awmean(x .± σx)
        @test wm.val ≈ 0.08737999999999996
        @test wm.err ≈ 0.07905694150420949
        @test m ≈ 38.44179426844445
        wm, m = gwmean(x .± σx)
        @test wm.val ≈ 0.08737999999999996
        @test wm.err ≈ 0.49016447665837415
        @test m ≈ 38.44179426844445
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

        # test chauvenet criterion
        x = [1.2, 1.5, 1.3, 2.4, 2.0, 2.1, 1.9, 2.2, 8.0, 2.3]
        xσ = [1,1,1,1,1,1,1,1,1,1.]
        expected = [true, true, true, true, true, true, true, true, false, true]
        @test Isoplot.chauvenet_func(x, xσ) == expected
        μ,σ,MSWD = wmean(x, xσ;chauvenet=true)
        @test μ ≈ 1.877777777777778
        @test MSWD ≈ 0.19444444444444445
        μ, MSWD = wmean((x .± xσ);chauvenet=true)
        @test value(μ) ≈ 1.877777777777778
        @test MSWD ≈ 0.19444444444444445
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

        wm, m = wmean(age68.(analyses[1:10]))
        @test wm.val ≈ 752.2453179272093
        @test wm.err ≈ 1.4781473739306696
        @test m ≈ 13.15644886325888

        m = mswd(age68.(analyses[1:10]))
        @test m ≈ 13.15644886325888
    end

    @testset "Concordia Metropolis" begin
        data = upperintercept.(0, analyses)
        @test Isoplot.dist_ll(ones(10), data, 751, 755) + Isoplot.prior_ll(data, 751, 755) ≈ -20.09136536048026
        @test Isoplot.dist_ll(ones(10), data, 750, 760) + Isoplot.prior_ll(data, 750, 760) ≈ -30.459633175497830
        @test Isoplot.dist_ll(ones(10), data, 752, 753) + Isoplot.prior_ll(data, 752, 753) ≈ -15.305463167234748
        @test Isoplot.dist_ll(ones(10), data, 751, 752) + Isoplot.prior_ll(data, 751, 752) ≈ -47.386667785224034

        tmindist, t0dist = metropolis_min(1000, ones(100), analyses; burnin=200, method=:bivariate)
        @test tmindist isa Vector{Float64}
        @test mean(tmindist) ≈ 751.85 atol = 1.5
        @test std(tmindist) ≈ 0.40 rtol = 0.6
        @test t0dist isa Vector{Float64}
        @test mean(t0dist) ≈ 80. atol = 90
        @test std(t0dist) ≈ 50. rtol = 0.6

        tmindist, tmaxdist, t0dist, lldist, acceptancedist = metropolis_minmax(10000, ones(10), analyses; burnin=200, method=:bivariate)
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

    # Try to import and calibrate a SIMS U-Pb dataset
    standards = importsimsdata("../examples/data")
    standardages = fill(1099., length(standards))
    calib = calibration(standards, standardages)

    @testset "Import" begin
        
        @test calib isa Isoplot.UPbSIMSCalibration{Float64}
        @test calib.line isa Isoplot.YorkFit{Float64}
        @test calib.line.xm ≈ 2.8160331226296806
        @test calib.line.ym.val ≈ 3.5574217244315136
        @test calib.line.ym.err ≈ 0.01166391028880385
        @test calib.line.slope.val ≈ 0.820053229510839
        @test calib.line.slope.err ≈ 0.0662908508944898
        @test calib.line.intercept.val ≈ 1.248124667809552
        @test calib.line.intercept.err ≈ 0.18704126735290513

        data = importsimsdata("../examples/data")
        @test age75(data[1], calib) ≈ [1113.4293385550338, 1053.0203346377668, 1034.4260533491415, 1054.7004885169204, 1068.990316049054, 1069.2081176470697, 1097.4490195779847, 1103.8607303163587, 1051.7684315423078, 1097.1247832545434, 1119.0003335172935, 1072.1944747704129, 1110.248360975909, 1127.7521530747485, 1087.7787393527146, 1158.1654711128456, 1069.4358908027727, 1048.2046253042013, 1125.243573300422, 1087.405486838316]
        @test age68(data[1], calib) ≈ [1101.7849182516, 1053.9899088381646, 1056.0073454530643, 1042.461945720766, 1019.4649813748944, 1061.9475493690124, 1083.06956873063, 1094.9459449280243, 1107.2169982835583, 1088.0301175750913, 1110.6013323577797, 1110.4392397549514, 1097.3607619368754, 1119.1362206959072, 1097.9650245711698, 1119.396230211389, 1111.2354657129576, 1070.4159730763058, 1093.654646933861, 1109.7764022094161] 

        analyses = calibrate(data, calib)
        @test analyses isa Vector{UPbAnalysis{Float64}}
        @test analyses[1] isa UPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === true
        @test value(age68(analyses[1])) ≈ 1086.114161562425
        @test stdev(age68(analyses[1])) ≈ 6.249296047616489
        @test isconcordant(analyses[10]) === false
        @test value(age68(analyses[10])) ≈ 1096.0020564376744
        @test stdev(age68(analyses[10])) ≈ 5.146380111453537
        @test value(age76(analyses[10])) ≈ 1089.7911572614648
        @test stdev(age76(analyses[10])) ≈ 10.373817886747076
        @test value(age75(analyses[10])) ≈ 1072.4641818116602
        @test stdev(age75(analyses[10])) ≈ 8.560315252726863
        @test value(ageconcordia(analyses[10])) ≈ 1093.6779441066528
        @test stdev(ageconcordia(analyses[10])) ≈ 5.040429228572438

        analyses = calibrate(data, calib, baseline=0.25)
        @test analyses isa Vector{UPbAnalysis{Float64}}
        @test analyses[1] isa UPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === true
        @test value(age68(analyses[1])) ≈ 1088.9985313345585
        @test stdev(age68(analyses[1])) ≈ 6.283919766775051
        @test isconcordant(analyses[10]) === true
        @test value(age68(analyses[10])) ≈ 1098.9097932521975
        @test stdev(age68(analyses[10])) ≈ 5.175351750349851
        @test value(age76(analyses[10])) ≈ 1097.5273085870058
        @test stdev(age76(analyses[10])) ≈ 4.395465527535868
        @test value(age75(analyses[10])) ≈ 1093.7737025147464
        @test stdev(age75(analyses[10])) ≈ 8.494053886674003
        @test value(ageconcordia(analyses[10])) ≈ 1098.4664884247334
        @test stdev(ageconcordia(analyses[10])) ≈ 5.09890800497881
    end

end

module PlotsTest

    using Test, Statistics
    using Measurements
    using Isoplot
    using ImageIO, FileIO, Plots

    import ..BaseTests: analyses, calib

    # Base.retry_load_extensions()
    @testset "Plots.jl Plotting" begin
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

        h = rankorder((1:10) .± (2*ones(10)), label="")
        savefig(h, "rankorder.png")
        img = load("rankorder.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9842757843137254,0.9882103758169932,0.9906218464052285) rtol = 0.01
        rm("rankorder.png")

        h = rankorder([CI(1:10) for _ in 1:10], label="")
        savefig(h, "CIa.png")
        img = load("CIa.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.973667565359477,0.9838965359477125,0.9901668300653594) rtol = 0.01
        rm("CIa.png")

        h = plot([CI(1:10) for _ in 1:10], label="")
        savefig(h, "CIb.png")
        img = load("CIb.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈  RGB{Float64}(0.9728901143790849,0.9831728594771241,0.9894754738562094) rtol = 0.01
        rm("CIb.png")

        h = plot([CI(1:10) for _ in 1:10], 1:10, label="")
        savefig(h, "CIc.png")
        img = load("CIc.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈  RGB{Float64}(0.9661086437908494,0.9805500816993463,0.9894015849673207) rtol = 0.01
        rm("CIc.png")

        h = plot([CI(1:10) for _ in 1:10], [CI(1:10) for _ in 1:10], label="")
        savefig(h, "CId.png")
        img = load("CId.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9829569444444444,0.9866916666666666,0.9890226307189541) rtol = 0.01
        rm("CId.png")

        # Plot SIMS calibration line
        h = plot(calib)
        savefig(h, "calib.png")
        img = load("calib.png")
        @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9293718627450989, 0.9393299346405233, 0.9497374999999998) rtol = 0.01
        rm("calib.png")
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
        @test sum(img)/length(img) ≈ RGBA{Float64}(0.8524913580246877, 0.8524913580246877, 0.9885884168482209, 1.0) rtol = 0.02
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
        @test sum(img)/length(img) ≈ RGBA{Float64}(0.9525158605664497, 0.9525158605664497, 0.9663422004357308, 1.0) rtol = 0.01
        rm("concordia.png")

        # Plot single concordia line
        f3 = Figure()
        ax3 = Axis(f3[1,1])
        concordialine!(0, 100)
        concordiacurve!(0,100)
        save("concordia.png",f3)
        img = load("concordia.png")
        @test size(img) == (900, 1200)
        @test sum(img)/length(img) ≈ RGBA{Float64}(0.9847984095860566, 0.9847984095860566, 0.9847984095860566, 1.0) rtol = 0.01
        rm("concordia.png")

    end
end

