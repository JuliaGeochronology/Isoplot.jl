module BaseTests

    using Test, Statistics, LinearAlgebra
    using Measurements, Distributions
    using Isoplot

    @testset "Show" begin
        yf = Isoplot.YorkFit(1±1, 1±1, 0.0, 1±1, 1.0)
        @test display(yf) != NaN

        ci = CI(1:10)
        @test ci == CI{Float64}(5.5, 3.0276503540974917, 5.5, 1.225, 9.775)
        @test "$ci" === "5.5 +4.28/-4.28"
        @test display(ci) != NaN

        a = Age(ci)
        @test display(a) != NaN

        mf = mcleanfit(1:10, ones(10), 1:10, ones(10))
        @test display(mf) != NaN
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

    @testset "U-Th-Pb" begin
        μ = [18.242, 0.59261, 0.16001]
        Σ = [0.00031 2.50e-6 2.50e-7; 
             2.50e-6 2.94e-8 2.50e-7;
             2.50e-7 2.50e-7 7.39e-7;]
        d1 = UThPbAnalysis(μ, Σ)
        @test d1 isa UThPbAnalysis{Float64}
        a75, a68, a82 = age(d1)
        @test a75.val ≈ 3000.0257585438353
        @test a75.err ≈ 1.9141367875881798
        @test a68.val ≈ 2999.985437721613
        @test a68.err ≈ 1.754454624116367
        @test a82.val ≈ 3000.0627649796165
        @test a82.err ≈ 21.190771017329745
        @test age68(d1) == a68
        @test age75(d1) == a75
        @test age82(d1) == a82

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

        # Generating an ellipse from an analysis
        e1 = Ellipse(d1, npoints=7)
        @test e1 isa Isoplot.Ellipse
        @test e1.x ≈ [22.70307499779583, 22.681635852743682, 22.63876085494785, 22.61732500220417, 22.638764147256317, 22.68163914505215, 22.70307499779583]
        @test e1.y ≈ [0.4089925091091875, 0.40901969166358015, 0.4086701825543926, 0.40829349089081246, 0.4082663083364198, 0.4086158174456074, 0.4089925091091875]

        # Generating an ellipse from means and correlation
        e2 = Ellipse(0., 1., 0., 1., 0., npoints=7)
        @test e2 isa Isoplot.Ellipse
        @test e2.x ≈ [0., -2.1198109374024336, -2.119810937402434, 0, 2.119810937402433, 2.1198109374024336, 0.]
        @test e2.y ≈ [2.447746830680816, 1.2238734153404085, -1.2238734153404074, -2.447746830680816, -1.2238734153404094, 1.223873415340408, 2.447746830680816]

        e3 = Ellipse(MvNormal([0.,0.], [1. 0.; 0. 1.]), npoints=7)
        @test e3 isa Isoplot.Ellipse
        @test e3.x ≈ [0., -2.1198109374024336, -2.119810937402434, 0, 2.119810937402433, 2.1198109374024336, 0.]
        @test e3.y ≈ [2.447746830680816, 1.2238734153404085, -1.2238734153404074, -2.447746830680816, -1.2238734153404094, 1.223873415340408, 2.447746830680816]

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
        @test stacey_kramers(0) == (18.7, 15.628, 38.63)
        @test stacey_kramers(1000) == (17.06558429108375, 15.509293514148855, 36.76149339195804)
        @test stacey_kramers(3700) == (11.152, 12.998, 31.23)
        @test stacey_kramers(4000) == (10.543938494456278, 12.30939102227753, 30.633643986839882)
        @test stacey_kramers(4567) == (9.314476625036953, 10.295372492782558, 29.482069098391527)
        @test stacey_kramers(5000) === (NaN, NaN, NaN)

    end
    @testset "Radiocarbon" begin

        d = Radiocarbon(1000, 10; calibration=intcal13)
        @test d isa Radiocarbon{Float64}
        @test d isa Distribution
        @test pdf(d, 950) ≈ 0.003534364129415577
        @test logpdf(d, 950) ≈ -5.645221874338272
        @test pdf.(d, [900, 950, 1000]) ≈ [3.458394021789428e-7, 0.003534364129415577, 1.4485272914446688e-12]
        @test cdf.(d, [-1, 900, 950, 1000, 1e5]) ≈ [0.0, 0.0030784606974826384, 0.9882404868727294, 0.9999999996129478, 1.0]
        @test ccdf.(d, [-1, 900, 950, 1000, 1e5]) ≈ [1.0, 0.9969218851419205, 0.01529387725668608, 3.8850078477011496e-10, 0.0]
        @test logcdf.(d, [-1, 900, 950, 1000, 1e5]) ≈ [-7649.088264850022, -5.7833255804696755, -0.011829203087212914, -3.8705216809375997e-10, -0.0]
        @test logccdf.(d, [-1, 900, 950, 1000, 1e5]) ≈ [-0.0, -0.003082861997615982, -4.180302709978416, -21.66872592612259, -8.706770436253262e7]

        @test eltype(d) === partype(d) === Float64
        @test location(d) ≈ 927.2556461305643
        @test scale(d) ≈ 7.5077650707247

        @test mean(d) ≈ 927.2556461305643
        @test var(d) ≈ 7.5077650707247^2
        @test std(d) ≈ 7.5077650707247
        @test skewness(d) ≈ -4.6327184250788696
        @test kurtosis(d) ≈ 67.68860469654331

        x = Radiocarbon(10000, 10, calibration=intcal13)
        @test mean(x) ≈ 11474.878002425487
        @test std(x) ≈ 82.19571126561112

        x = Radiocarbon(10000, 10)
        @test mean(x) ≈ 11476.343648266482
        @test std(x) ≈ 86.70372094312096

        age_14C = [9962.71, 10002.18, 10007.3, 10012.88, 10016.1, 10019.81, 10026.88, 10047.97, 10056.09, 10065.81]
        age_sigma_14C = fill(25, 10)
        ages = Radiocarbon.(age_14C, age_sigma_14C, calibration=intcal20)
        @test ages isa Vector{<:Radiocarbon}

        # Weighted means on radiocarbon ages
        wm, mswd = wmean(ages)
        @test value(wm) ≈ 11509.329228181394
        @test stdev(wm) ≈ 30.825739204181755
        @test mswd ≈ 0.3705296426717654
        
        # Monte Carlo weighted means on distributions
        distwm = distwmean(ages, N=50000)
        @test mean(distwm) ≈ 11509.2 atol = 1
        @test std(distwm) ≈ 47.4 atol = 1

        # Metropolis on radiocarbon ages
        tmindist = metropolis_min(5000, ones(10), ages; burnin=200)
        @test mean(tmindist) ≈ 11453 atol=25
        @test std(tmindist) ≈ 67.8 atol=35

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
        @test value(c.line.slope) ≈ slope           atol=4stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.706       atol=4stdev(c.line.intercept)
        @test value(age(c)) ≈ 3502.5243481286466    atol=4stdev(age(c))
        @test stdev(age(c)) ≈ 5.104097410724326     atol=1
        @test value(age(slope, Isoplot.λ87Rb)) ≈ 3502.5243481286466

        analyses = [LuHfAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.2828) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, LuHfAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=4stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.2828      atol=4stdev(c.line.intercept)
        @test value(age(c)) ≈ 2613.2921354810956    atol=4stdev(age(c))
        @test stdev(age(c)) ≈ 5.6338908837666954    atol=2
        @test value(age(slope, Isoplot.λ176Lu)) ≈ 2613.2921354810956

        analyses = [ReOsAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.14848) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, ReOsAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=4stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.14848     atol=4stdev(c.line.intercept)
        @test value(age(c)) ≈ 2923.492370389601     atol=4stdev(age(c))
        @test stdev(age(c)) ≈ 2.7958914719906987    atol=1
        @test value(age(slope, Isoplot.λ187Re)) ≈ 2923.492370389601

        analyses = [SmNdAnalysis(0.01randn(N).+x, 0.001randn(N).+slope.*x.+0.5126) for x in 1:10]
        c = isochron(analyses)
        @test c isa Isoplot.Isochron{Float64, SmNdAnalysis{Float64}}
        @test value(c.line.slope) ≈ slope           atol=4stdev(c.line.slope)
        @test value(c.line.intercept) ≈ 0.5126      atol=4stdev(c.line.intercept)
        @test value(age(c)) ≈ 7478.565936454944     atol=4stdev(age(c))
        @test stdev(age(c)) ≈ 13.869217589658614    atol=3
        @test value(age(slope, Isoplot.λ147Sm)) ≈ 7478.565936454944
    end

    @testset "Weighted means" begin
        # Using vectors of mean and standard deviation
        x = [-3.4699, -0.875, -1.4189, 1.2993, 1.1167, 0.8357, 0.9985, 1.2789, 0.5446, 0.5639]
        σx = ones(10)/4
        @test all(awmean(x, σx) .≈ (0.08737999999999996, 0.07905694150420949, 38.44179426844445))
        @test all(gwmean(x, σx) .≈ (0.08737999999999996, 0.49016447665837415, 38.44179426844445))
        # Using Measurements
        wm, m = wmean(x .± σx; corrected=false)
        @test value(wm) ≈ 0.08737999999999996
        @test stdev(wm) ≈ 0.07905694150420949
        @test m ≈ 38.44179426844445
        wm, m = wmean(x .± σx; corrected=true)
        @test value(wm) ≈ 0.08737999999999996
        @test stdev(wm) ≈ 0.49016447665837415
        @test m ≈ 38.44179426844445
        # Using Distributions
        wm, m = wmean(Normal.(x, σx); corrected=false)
        @test value(wm) ≈ 0.08737999999999996
        @test stdev(wm) ≈ 0.07905694150420949
        @test m ≈ 38.44179426844445
        wm, m = wmean(Normal.(x, σx); corrected=true)
        @test value(wm) ≈ 0.08737999999999996
        @test stdev(wm) ≈ 0.49016447665837415
        @test m ≈ 38.44179426844445
        # MSWDs
        @test mswd(x, σx) ≈ 38.44179426844445
        @test mswd(x .± σx) ≈ 38.44179426844445
        @test mswd(Normal.(x, σx)) ≈ 38.44179426844445

        # Test underdispersed data
        σx .*= 20 
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
        @test mean(d) ≈ μ atol = 0.3
        @test std(d) ≈ σ atol = 3

        d = distwmean(a,b,c; corrected=false)
        μ,σ,_ = wmean([0,50,100], [2,1,3]; corrected=false)
        @test mean(d) ≈ μ atol = 0.2
        @test std(d) ≈ σ atol = 0.02
        d = distwmean(a,b,c; corrected=true)
        μ,σ,_ = wmean([0,50,100], [2,1,3]; corrected=true)
        @test mean(d) ≈ μ atol = 0.3
        @test std(d) ≈ σ atol = 3

        # test Chauvenet's criterion
        x = [1.2, 1.5, 1.3, 2.4, 2.0, 2.1, 1.9, 2.2, 8.0, 2.3]
        xσ = [1,1,1,1,1,1,1,1,1,1.]
        expected = [true, true, true, true, true, true, true, true, false, true]
        @test (Isoplot.chauvenet_criterion(x, xσ) .> 1) == expected
        μ,σ,MSWD = wmean(chauvenet(x, xσ)...)
        @test μ ≈ 1.877777777777778
        @test MSWD ≈ 0.19444444444444445
        μ, MSWD = wmean(chauvenet(x .± xσ))
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
        # CIs
        ci = CI(randn(10_000))
        @test ci.mean ≈ 0 atol = 0.1
        @test ci.sigma ≈ 1 atol = 0.1
        @test ci.lower ≈ -1.9599639845400576 atol = 0.15
        @test ci.upper ≈ 1.9599639845400576 atol = 0.15
        
        ci = CI([0.4227751675791465, 2.339746777712491, 0.7073822332289744, 1.115676532764763, 1.3295183295746795, 1.1839524860728812, 9.051763176763613, 3.163690972571135, 0.35529944660568286, 20.552275066037705, 1.3021911627717235, 0.8154804384732102, 2.630294225259665, 0.24662095224118538, 0.4755015782877633])
        @test ci ≈ CI{Float64}(3.0461445697296416, 5.314497528659032,1.1839524860728812, 0.2846584252687595, 16.52709590479178)

        # Simple linear regression
        ϕ = lsqfit(1:10, 1:10)
        @test ϕ[1] ≈ 0 atol=1e-12
        @test ϕ[2] ≈ 1 atol=1e-12

        # Tera-Wasserburg conversion
        d = UPbAnalysis(0.6, 0.006, 0.08, 0.0008, 0.7)
        t = twanalysis(d)
        @test t isa Analysis2D{Float64}
        @test t.μ ≈ [1/0.08, 0.6/0.08/value(Isoplot.U238_U235)]
        samples = rand(d, 1_000_000) # Compare to Monte Carlo
        X = 1 ./ samples[2,:]
        Y = samples[1,:] ./ samples[2,:] ./ value(Isoplot.U238_U235)
        @test t.σ ≈ [std(X), std(Y)] rtol=0.02
        @test t.Σ[1,2]/prod(t.σ) ≈ cor(X, Y) atol=0.02
        @test t.Σ == t.Σ'
        @test twanalysis(UThPbAnalysis([0.6, 0.08, 0.02], [0.006, 0.0008, 0.0002], [0.006^2 0.7*0.006*0.0008 0; 0.7*0.006*0.0008 0.0008^2 0; 0 0 0.0002^2])).μ ≈ t.μ
        # Concordia in Tera-Wasserburg space
        x, y = Isoplot.concordiacoords(1000, value(λ235U), value(λ238U), tw=true)
        @test x ≈ 1/Isoplot.ratio(1000, value(λ238U))
        @test y ≈ Isoplot.ratio(1000, value(λ235U))/Isoplot.ratio(1000, value(λ238U))/value(Isoplot.U238_U235)
        @test Isoplot.ageticks((0, 1000)) == [0.0, 500.0, 1000.0]

        # York (1968) fit
        x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
        y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
        σx = σy = ones(10)/4
        yf = yorkfit(x, σx, y, σy)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈-0.2446693790977319
        @test yf.intercept.err ≈ 0.2469541320914601
        @test yf.slope.val ≈ 1.0428730084538775
        @test yf.slope.err ≈ 0.039561084436542084
        @test yf.mswd ≈ 1.4272439423337884

        yf = yorkfit(x, σx, y, σy, zeros(length(x)))
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈-0.2446693790977319
        @test yf.intercept.err ≈ 0.2469541320914601
        @test yf.slope.val ≈ 1.0428730084538775
        @test yf.slope.err ≈ 0.039561084436542084
        @test yf.mswd ≈ 1.4272439423337884

        # `corrected` scales uncertainties by √max(MSWD,1), and nothing else
        yc = yorkfit(x, σx, y, σy, zeros(length(x)); corrected=true)
        yu = yorkfit(x, σx, y, σy, zeros(length(x)))
        f = sqrt(max(yu.mswd, 1))
        @test yc.slope.val == yu.slope.val
        @test yc.mswd == yu.mswd
        @test yc.slope.err ≈ f * yu.slope.err
        @test yc.intercept.err ≈ f * yu.intercept.err
        @test yc.ym.err ≈ f * yu.ym.err
        # Underdispersed data are not shrunk
        ys = yorkfit(x, 10σx, y, 10σy, zeros(length(x)); corrected=true)
        @test ys.slope.err == yorkfit(x, 10σx, y, 10σy, zeros(length(x))).slope.err
        # Agrees with mcleanfit
        mc = mcleanfit(x, σx, y, σy; corrected=true)
        @test yc.slope.err ≈ mc.slope.err rtol=0.01

        x = ((1:100) .+ randn.()) .± 1
        y = (2*(1:100) .+ randn.()) .± 1
        yf = yorkfit(x, y)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈ 0 atol = 2
        @test yf.slope.val ≈ 2 atol = 0.1
        @test yf.mswd ≈ 1 atol = 0.5

        yf = yorkfit(analyses)
        @test yf isa Isoplot.YorkFit
        @test yf.intercept.val ≈ 0.004516874477465516 
        @test yf.intercept.err ≈ 0.0030191000542899014
        @test yf.slope.val ≈ 0.10849669738014765
        @test yf.slope.err ≈ 0.0027553160759274534
        @test yf.xm ≈ 1.095717612358936
        @test yf.ym.val ≈ 0.12339861666825475
        @test yf.ym.err ≈ 1.7668787238429157e-5
        @test yf.mswd ≈ 0.665079939690809

        ui = upperintercept(analyses)
        @test ui.val ≈ 752.5343109308996 
        @test ui.err ≈ 0.4124832290832207

        li = lowerintercept(analyses)
        @test li.val ≈ 103.1132456354392
        @test li.err ≈ 69.42485829479817

        ui, li = intercepts(analyses)
        @test ui.val ≈ 752.5343109308996
        @test ui.err ≈ 0.4124832290832207
        @test li.val ≈ 103.1132456354392
        @test li.err ≈ 69.42485829479817

        wm, m = wmean(age68.(analyses[1:10]))
        @test wm.val ≈ 752.2453179272093
        @test wm.err ≈ 1.4781473739306696
        @test m ≈ 13.15644886325888

        m = mswd(age68.(analyses[1:10]))
        @test m ≈ 13.15644886325888
    end

    @testset "McLean fit" begin
        x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
        y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
        σx = σy = ones(10)/4
        r = [0.0, 0.3, -0.4, 0.0, 0.5, 0.1, 0.0, -0.2, 0.2, 0.0]

        # The fitted line should be the same as for the York fit, in 2D
        for ρ in (zeros(10), r)
            mf, yf = mcleanfit(x, σx, y, σy, ρ), yorkfit(x, σx, y, σy, ρ)
            @test mf isa Isoplot.McLeanFit{Float64}
            @test mf.intercept.val ≈ yf.intercept.val rtol=1e-6
            @test mf.slope.val ≈ yf.slope.val rtol=1e-6
            @test mf.mswd ≈ yf.mswd rtol=1e-6
            # The points where height and slope are uncorrelated can differ, but should be on the same line
            @test Isoplot.line(mf, yf.xm).val ≈ yf.ym.val rtol=1e-6
            # Uncertainties are from the exact Hessian rather than York's approximation
            @test mf.intercept.err ≈ yf.intercept.err rtol=1e-2
            @test mf.slope.err ≈ yf.slope.err rtol=1e-2
            @test mf.ym.err ≈ yf.ym.err rtol=1e-2
        end
        mf = mcleanfit(x, σx, y, σy)
        @test mf.intercept.val ≈ -0.2446693790977319 rtol=1e-6
        @test mf.slope.val ≈ 1.0428730084538775 rtol=1e-6
        @test mf.mswd ≈ 1.4272439423337884 rtol=1e-6

        # Line is the right height at xm, where uncertainties of height and slope are uncorrelated
        @test Isoplot.line(mf, mf.xm) == mf.ym
        @test Isoplot.line(mf, 0.0).val ≈ mf.intercept.val atol=1e-8
        @test Isoplot.invline(mf, Isoplot.line(mf, 3.0).val) ≈ 3.0

        # Expanded uncertainties
        mfc = mcleanfit(x, σx, y, σy, corrected=true)
        @test mfc.slope.val == mf.slope.val
        @test mfc.slope.err ≈ mf.slope.err * sqrt(mf.mswd)
        @test mfc.intercept.err ≈ mf.intercept.err * sqrt(mf.mswd)
        @test mcleanfit(x, σx, 10σy, 10σy, corrected=true).slope.err ≈ mcleanfit(x, σx, 10σy, 10σy).slope.err # mswd < 1

        # Analytic gradient and Hessian against finite differences, in two and three dimensions
        for N in (2, 3)
            p = [Isoplot.SVector{N}(ntuple(k -> (i, 2i + sin(i), -i + cos(3i))[k], N)) for i in 1:8]
            Σ = [(A = Isoplot.SMatrix{N,N}(cos.((1:N^2) .+ i)); A*A' + 0.1*Isoplot.SMatrix{N,N}(I)) for i in 1:8]
            W = inv.(Σ)
            z, h = Isoplot.SVector{2N}(ntuple(k -> (1.0, 0.5, 2.0, -0.3, 1.0, 0.8)[k], 2N)), 1e-5
            L(z) = Isoplot._mclean_loglikelihood(Isoplot._mclean_halves(z, Val(N))..., p, W)
            G(z) = first(Isoplot._mclean_derivatives(Isoplot._mclean_halves(z, Val(N))..., p, W))
            e(k) = Isoplot.SVector{2N}(ntuple(i -> i == k ? h : 0.0, 2N))
            Gz, Hz = Isoplot._mclean_derivatives(Isoplot._mclean_halves(z, Val(N))..., p, W)
            for k in 1:2N
                @test Gz[k] ≈ (L(z + e(k)) - L(z - e(k)))/2h rtol=1e-6
                @test Hz[:,k] ≈ (G(z + e(k)) - G(z - e(k)))/2h rtol=1e-6
            end
        end

        # Other input types
        @test mcleanfit(Tuple(x), Tuple(σx), Tuple(y), Tuple(σy)).slope.val ≈ mf.slope.val
        @test mcleanfit(Float32.(x), Float32.(σx), Float32.(y), Float32.(σy)) isa Isoplot.McLeanFit{Float32}
        @test mcleanfit(Float32.(x), Float32.(σx), Float32.(y), Float32.(σy)).slope.val ≈ mf.slope.val rtol=1e-3
        @test mcleanfit(x .± 0.25, y .± 0.25).slope.val ≈ mf.slope.val rtol=1e-6
        mf = mcleanfit(analyses)
        yf = yorkfit(analyses)
        @test mf.slope.val ≈ yf.slope.val rtol=1e-6
        @test mf.intercept.val ≈ yf.intercept.val rtol=1e-4
        @test Isoplot.line(mf, yf.xm).val ≈ yf.ym.val rtol=1e-6
        @test mf.mswd ≈ yf.mswd rtol=1e-6

        # NaNs are ignored
        xn = copy(x); xn[3] = NaN
        keep = [1:2; 4:10]
        @test mcleanfit(xn, σx, y, σy).slope.val ≈ mcleanfit(x[keep], σx[keep], y[keep], σy[keep]).slope.val

        # Perfect fit, which has zero chi-squared and so a log-likelihood of zero
        mf = mcleanfit(1.0:6.0, fill(0.1, 6), 2.0 .* (1:6) .+ 1, fill(0.1, 6))
        @test mf.slope.val ≈ 2
        @test mf.intercept.val ≈ 1
        @test mf.mswd ≈ 0 atol=1e-10

        # Fits in three dimensions
        N, n = 3, 20
        a₀, v₀ = [0.5, 1.0, 2.0], [1.0, 0.5, -2.0]
        covs = [(A = [1 0.3 -0.2; 0.1 1 0.4; 0.2 -0.3 1] * (0.05 + 0.1*i/n); A*A') for i in 1:n]
        ts = range(1, 10, length=n)
        mf = mcleanfit([a₀ .+ v₀*t for t in ts], covs) # Exact line
        @test mf isa Isoplot.McLeanFit{Float64,3}
        @test mf.v ≈ v₀
        @test mf.mswd ≈ 0 atol=1e-10
        @test mf.slope isa Isoplot.SVector{2,<:Measurement}
        @test mf.intercept[1].val ≈ a₀[2] - a₀[1]*v₀[2]
        @test mf.intercept[2].val ≈ a₀[3] - a₀[1]*v₀[3]
        @test_throws ArgumentError mf.ym

        # Statistical behavior given noisy data: slope z-scores and MSWD should be ~standard
        zs, mswds = Float64[], Float64[]
        for _ in 1:200
            pts = [a₀ .+ v₀*ts[i] .+ cholesky(Symmetric(covs[i])).L * randn(3) for i in 1:n]
            mf = mcleanfit(pts, covs)
            append!(zs, (mf.v[2:3] .- v₀[2:3]) ./ sqrt.(diag(mf.Σ)[3:4]))
            push!(mswds, mf.mswd)
        end
        @test std(zs) ≈ 1 atol=0.2
        @test abs(mean(zs)) < 0.2
        @test mean(mswds) ≈ 1 atol=0.1

        # Same fit regardless of the order of the other variables
        pts = [a₀ .+ v₀*ts[i] .+ cholesky(Symmetric(covs[i])).L * randn(3) for i in 1:n]
        # Variables of the fitted line, with uncertainty, at a given value of the first variable
        mf2 = mcleanfit(1.0:6.0, fill(0.1, 6), 2.0 .* (1:6) .+ 1, fill(0.1, 6))
        for x in (0.0, 3.5, 10.0)
            @test Isoplot.line(mf2, x, 2).val ≈ Isoplot.line(mf2, x).val
            @test Isoplot.line(mf2, x, 2).err ≈ Isoplot.line(mf2, x).err
        end
        @test_throws ArgumentError Isoplot.line(mf2, 1.0, 3)
        yf = yorkfit(1.0:6.0, fill(0.1, 6), 2.0 .* (1:6) .+ 1, fill(0.1, 6))
        @test Isoplot.line(yf, 2.0, 2) == Isoplot.line(yf, 2.0)
        @test_throws ArgumentError Isoplot.line(yf, 2.0, 3)
        mf3 = mcleanfit(pts, covs)
        # Ellipse of a pair of variables, and ellipsoid, of a 3D analysis
        a3 = Analysis([1.0, 2.0, 3.0], [0.2, 0.3, 0.4], [0.04 0.03 0.0; 0.03 0.09 0.05; 0.0 0.05 0.16])
        e13 = Ellipse(a3, 1, 3)
        @test (e13.x₀, e13.y₀, e13.σx₀, e13.σy₀) == (1.0, 3.0, 0.2, 0.4)
        @test e13.x == Ellipse(Analysis2D([1.0, 3.0], [0.2, 0.4], [0.04 0.0; 0.0 0.16])).x
        x, y, z = Isoplot.ellipsoid(a3, sigmalevel=2, npoints=10)
        @test size(x) == size(y) == size(z) == (10, 10)
        Σ⁻¹ = inv(Matrix(a3.Σ))
        @test all(isapprox(([x[i]-1, y[i]-2, z[i]-3]' * Σ⁻¹ * [x[i]-1, y[i]-2, z[i]-3]), 4) for i in eachindex(x))
        a, v = mf3.a, mf3.v
        @test Isoplot.line(mf3, 1.5, 3).val ≈ a[3] + v[3]*(1.5 - a[1])
        @test Isoplot.line(mf3, a[1], 2).err ≈ sqrt(mf3.Σ[1,1]) # At a₁ only the uncertainty of a₂ contributes
        @test Isoplot.line(mf3, 1.5, 3).err > 0
        @test_throws ArgumentError Isoplot.line(mf3, 1.5, 1)
        @test_throws ArgumentError Isoplot.line(mf3, 1.5, 4)

        mf = mcleanfit(pts, covs)
        mfp = mcleanfit([q[[1,3,2]] for q in pts], [C[[1,3,2],[1,3,2]] for C in covs])
        @test mfp.v ≈ mf.v[[1,3,2]]
        @test mfp.a ≈ mf.a[[1,3,2]]
        @test mfp.mswd ≈ mf.mswd
        @test mfp.Σ ≈ mf.Σ[[2,1,4,3],[2,1,4,3]]
        @test mcleanfit(pts, covs, corrected=true).Σ ≈ mf.Σ * max(mf.mswd, 1)

        # From `Analysis` objects, and with NaNs
        @test mcleanfit([Analysis(pts[i], covs[i]) for i in 1:n]).v ≈ mf.v
        pts[4] = [NaN, 1, 1]
        @test mcleanfit(pts, covs).v ≈ mcleanfit(pts[[1:3; 5:n]], covs[[1:3; 5:n]]).v

        # Robustness with noisy data, where the fit is poorly determined: never fails to converge or
        # settles on a saddle point, and is at least as good as the York fit (which can fail to converge)
        for _ in 1:33
            n = 12
            σx, σy, r = 0.1 .+ 2rand(n), 0.1 .+ 2rand(n), 1.6rand(n) .- 0.8
            x = 10rand(n)
            y = 0.7 .* x .+ 1 .+ randn(n)
            mf = mcleanfit(x, σx, y, σy, r)
            @test !isnan(mf.mswd)
            W = [inv(Isoplot.SMatrix{2,2}(σx[i]^2, r[i]*σx[i]*σy[i], r[i]*σx[i]*σy[i], σy[i]^2)) for i in 1:n]
            p = [Isoplot.SVector(x[i], y[i]) for i in 1:n]
            L(f) = Isoplot._mclean_loglikelihood(Isoplot.SVector(0.0, f.intercept.val), Isoplot.SVector(1.0, f.slope.val), p, W)
            @test L(mf) >= L(yorkfit(x, σx, y, σy, r, maxiter=100)) - 1e-9
        end
        for _ in 1:33
            covs = [(A = randn(3,3)/2 + I; 0.03*A*A') for _ in 1:15]
            @test !isnan(mcleanfit([10rand()*[1, -0.3, 0.6] + cholesky(Symmetric(C)).L*randn(3) for C in covs], covs).mswd)
        end

        # Degenerate cases
        @test isnan(mcleanfit([1.0], [1.0], [1.0], [1.0]).slope)
        @test isnan(mcleanfit([1.0, 2.0], [1.0, 1.0], [1.0, 2.5], [1.0, 1.0]).mswd) # no degrees of freedom
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
        @test Isoplot.dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100, 101) ≈ -3.693335779335658
        @test Isoplot.dist_ll(MeltsVolcanicZirconDistribution, Normal.(mu, sigma), 100, 101) ≈ -4.736451831031247

        tmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03

        tmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu .± sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03

        tmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, Normal.(mu, sigma), burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03

        tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03
        @test mean(tmaxdist) ≈ 101.08  atol=0.015
        @test std(tmaxdist) ≈ 0.06917  atol=0.03
        @test mean(acceptancedist) ≈ 0.6 atol=0.2
        @test acceptancedist isa BitVector
        @test lldist isa Vector{Float64}

        tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu .± sigma, burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03
        @test mean(tmaxdist) ≈ 101.08  atol=0.015
        @test std(tmaxdist) ≈ 0.06917  atol=0.03
        @test mean(acceptancedist) ≈ 0.6 atol=0.2

        tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, Normal.(mu, sigma), burnin=10^5)
        @test mean(tmindist) ≈ 99.9228  atol=0.015
        @test std(tmindist) ≈ 0.11355  atol=0.03
        @test mean(tmaxdist) ≈ 101.08  atol=0.015
        @test std(tmaxdist) ≈ 0.06917  atol=0.03
        @test mean(acceptancedist) ≈ 0.6 atol=0.2

        @test mean(UniformDistribution) ≈ 1
        @test mean(TriangularDistribution) ≈ 1
        @test mean(HalfNormalDistribution) ≈ 1
        @test mean(ExponentialDistribution) ≈ 1.03 atol=0.01
        @test mean(MeltsZirconDistribution) ≈ 1 atol=0.01
        @test mean(MeltsVolcanicZirconDistribution) ≈ 1 atol=0.01
    end

    @testset "Import" begin
        # Try to import and calibrate a SIMS U-Pb dataset
        standards = importsimsdata("../examples/data")
        standardages = fill(1099., length(standards))

        calib = calibration(standards, standardages, correction=:Pb204)
        @test calib isa Isoplot.UPbSIMSCalibration{Float64}
        @test calib.line isa Isoplot.YorkFit{Float64}
        @test calib.line.xm ≈ 2.8431334173416767
        @test calib.line.ym.val ≈ 3.5857266263911396
        @test calib.line.ym.err ≈ 0.0066776603378822415
        @test calib.line.slope.val ≈ 0.8050681793138171
        @test calib.line.slope.err ≈ 0.03899643939328304
        @test calib.line.intercept.val ≈ 1.2968103825456048
        @test calib.line.intercept.err ≈ 0.11107299073272647
        @test calib.line.mswd ≈ 0.9312416907930794 

        calib = calibration(standards, standardages, correction=:Pb208)
        @test calib isa Isoplot.UPbSIMSCalibration{Float64}
        @test calib.line isa Isoplot.YorkFit{Float64}
        @test calib.line.xm ≈ 2.763415189906566
        @test calib.line.ym.val ≈ 3.5802230412292295
        @test calib.line.ym.err ≈ 0.006750554163821435
        @test calib.line.slope.val ≈ 0.7871464039573363
        @test calib.line.slope.err ≈ 0.03912661929741835
        @test calib.line.intercept.val ≈ 1.4050107118531958
        @test calib.line.intercept.err ≈ 0.1083336210898975
        @test calib.line.mswd ≈ 1.428831456701787

        calib = calibration(standards, standardages, correction=:none)
        @test calib isa Isoplot.UPbSIMSCalibration{Float64}
        @test calib.line isa Isoplot.YorkFit{Float64}
        @test calib.line.xm ≈ 2.8461846235947252 
        @test calib.line.ym.val ≈ 3.585253321652293
        @test calib.line.ym.err ≈ 0.006691503076240842
        @test calib.line.slope.val ≈ 0.8098407435571999
        @test calib.line.slope.err ≈ 0.039102112286309
        @test calib.line.intercept.val ≈ 1.2802970497792723
        @test calib.line.intercept.err ≈ 0.11149281502742348
        @test calib.line.mswd ≈ 0.8789962702237475 

        # Test import with baseline subtraction
        data = importsimsdata("../examples/data", baseline=0.25)
        analyses = calibrate(data, calib)
        @test analyses isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === true
        @test value(age68(analyses[1])) ≈ 1089.1405666597918
        @test stdev(age68(analyses[1])) ≈ 6.273971991547559
        @test isconcordant(analyses[10]) === true
        @test value(age68(analyses[10])) ≈ 1101.2698237147065
        @test stdev(age68(analyses[10])) ≈ 5.18555607544457
        @test value(age76(analyses[10])) ≈ 1099.669945482014 
        @test stdev(age76(analyses[10])) ≈ 4.402729107204171
        @test value(age75(analyses[10])) ≈ 1095.3362790870392
        @test stdev(age75(analyses[10])) ≈ 8.500823136680218
        @test value(ageconcordia(analyses[10])) ≈ 1100.752826666941
        @test stdev(ageconcordia(analyses[10])) ≈ 5.1073391063212235
        @test value(ageconcordia3d(analyses[10])) ≈ 1096.507620343022
        @test stdev(ageconcordia3d(analyses[10])) ≈ 5.046405074810369

        # Test import without baseline subtraction
        data = importsimsdata("../examples/data")
        @test age75(data[1], calib) ≈ [1114.5967327821122, 1054.1518772448574, 1035.546123711327, 1055.8330574133802, 1070.131545207158, 1070.3494778604183, 1098.607136588629, 1105.0225871551654, 1052.899208342237, 1098.2827105154104, 1120.1709277117836, 1073.337629093534, 1111.4139201551982, 1128.9277389239833, 1088.931170951752, 1159.3580721712417, 1070.5773880413249, 1049.3332168881145, 1126.4177327678221, 1088.5576979036687]
        @test age68(data[1], calib) ≈ [1103.5347418354286, 1055.6698782743545, 1057.6902739428133, 1044.1249888384261, 1021.094167739909, 1063.6391852056834, 1084.7921009177614, 1096.6858048734648, 1108.9747282758387, 1089.7598911435075, 1112.363984885957, 1112.201656577389, 1099.104141216443, 1120.9112757873825, 1099.7092842915501, 1121.1716628803867, 1112.999040304355, 1072.12000832304, 1095.3926244105648, 1111.5378551098597]

        analyses = calibrate(data, calib)
        @test analyses isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === true
        @test value(age68(analyses[1])) ≈ 1087.841138890462
        @test stdev(age68(analyses[1])) ≈ 6.258425035757806
        @test isconcordant(analyses[10]) === false
        @test value(age68(analyses[10])) ≈ 1098.3577219643148
        @test stdev(age68(analyses[10])) ≈ 5.156537812728746
        @test value(age76(analyses[10])) ≈ 1091.9223110317737
        @test stdev(age76(analyses[10])) ≈ 10.736362149559248
        @test value(age75(analyses[10])) ≈ 1074.0105609640939
        @test stdev(age75(analyses[10])) ≈ 8.567283982507155
        @test value(ageconcordia(analyses[10])) ≈ 1095.9294600191558
        @test stdev(ageconcordia(analyses[10])) ≈ 5.046885422306974
        @test value(ageconcordia3d(analyses[10])) ≈ 1093.1104415499756
        @test stdev(ageconcordia3d(analyses[10])) ≈ 4.986436132488419

        analyses = calibrate(data, calib, correction=:Pb208)
        @test analyses isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === false
        @test value(age68(analyses[1])) ≈ 1053.7229717728435
        @test stdev(age68(analyses[1])) ≈ 6.261712417143128
        @test isconcordant(analyses[10]) === false
        @test value(age68(analyses[10])) ≈ 1066.9938629249943
        @test stdev(age68(analyses[10])) ≈ 5.225172489705245
        @test value(age76(analyses[10])) ≈ 1030.3032786714866
        @test stdev(age76(analyses[10])) ≈ 90.13006908470071
        @test value(age75(analyses[10])) ≈ 808.8996265223315
        @test stdev(age75(analyses[10])) ≈ 12.76612891824328
        @test value(ageconcordia(analyses[10])) ≈ 1060.933644369098 
        @test stdev(ageconcordia(analyses[10])) ≈ 5.026431992041755
        @test isnan(value(ageconcordia3d(analyses[10])))
        @test isnan(stdev(ageconcordia3d(analyses[10])))

        analyses = calibrate(data, calib, correction=:none)
        @test analyses isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1]) === true
        @test value(age68(analyses[1])) ≈ 1089.231227194806
        @test stdev(age68(analyses[1])) ≈ 6.2582952444563364
        @test isconcordant(analyses[10]) === true
        @test value(age68(analyses[10])) ≈ 1101.7668027149605
        @test stdev(age68(analyses[10])) ≈ 5.09389330557594
        @test value(age76(analyses[10])) ≈ 1100.9219832103886
        @test stdev(age76(analyses[10])) ≈ 4.173013425616315
        @test value(age75(analyses[10])) ≈ 1099.1242550059007
        @test stdev(age75(analyses[10])) ≈ 7.405826584618455
        @test value(ageconcordia(analyses[10])) ≈ 1101.2742458636797
        @test stdev(ageconcordia(analyses[10])) ≈ 4.906484873056218
        @test value(ageconcordia3d(analyses[10])) ≈ 1122.4637364497573
        @test stdev(ageconcordia3d(analyses[10])) ≈ 4.8200469677696

        analyses = calibrate_blockwise(data, calib, blocksize=5, correction=:Pb204)
        @test analyses isa Vector{Vector{UThPbAnalysis{Float64}}}
        @test analyses[1] isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1][1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1][1]) === true
        @test value(age68(analyses[1][1])) ≈ 1134.9098618894166
        @test stdev(age68(analyses[1][1])) ≈ 14.456203544532643
        @test isconcordant(analyses[10][1]) === false
        @test value(age68(analyses[10][1])) ≈ 1168.6844678115342
        @test stdev(age68(analyses[10][1])) ≈ 5.504561837703628
        @test value(age76(analyses[10][1])) ≈ 1165.2978111991247
        @test stdev(age76(analyses[10][1])) ≈ 13.882539521778183
        @test value(age75(analyses[10][1])) ≈ 1108.390682537911
        @test stdev(age75(analyses[10][1])) ≈ 22.425563301735476
        @test value(ageconcordia(analyses[10][1])) ≈ 1175.1301425017025
        @test stdev(ageconcordia(analyses[10][1])) ≈ 5.061803558541107
        @test value(ageconcordia3d(analyses[1][1])) ≈ 1194.3710016032248
        @test stdev(ageconcordia3d(analyses[1][1])) ≈ 12.154638963214998

        analyses = calibrate_blockwise(data, calib, blocksize=5, correction=:Pb208)
        @test analyses isa Vector{Vector{UThPbAnalysis{Float64}}}
        @test analyses[1] isa Vector{UThPbAnalysis{Float64}}
        @test analyses[1][1] isa UThPbAnalysis{Float64}
        @test isconcordant(analyses[1][1]) === false
        @test value(age68(analyses[1][1])) ≈ 1098.1481738040404
        @test stdev(age68(analyses[1][1])) ≈ 15.391305095639025
        @test isconcordant(analyses[10][1]) === false
        @test value(age68(analyses[10][1])) ≈ 1134.4105650817942
        @test stdev(age68(analyses[10][1])) ≈ 4.469372331888714
        @test value(age76(analyses[10][1])) ≈ 1112.5632113093047
        @test stdev(age76(analyses[10][1])) ≈ 79.59877929401242
        @test value(age75(analyses[10][1])) ≈ 822.5525431488766
        @test stdev(age75(analyses[10][1])) ≈ 16.13823147769109
        @test value(ageconcordia(analyses[10][1])) ≈ 1114.8534888779454
        @test stdev(ageconcordia(analyses[10][1])) ≈ 4.2928903087922405
        @test isnan(value(ageconcordia3d(analyses[1][1])))
        @test isnan(stdev(ageconcordia3d(analyses[1][1])))
    end

    # Import and calibrate a SIMS U-Pb dataset, for later use
    standards = importsimsdata("../examples/data")
    standardages = fill(1099., length(standards))
    calib = calibration(standards, standardages, correction=:none)
end

module PlotsTest

    using Test, Statistics
    using LinearAlgebra: Diagonal
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
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.8151617156862782,0.8151617156862782,0.986395212418301) rtol = 0.02
        rm("concordia.png")

        # Plot many concordia ellipses and concordia curve
        h = plot(analyses, color=:blue, alpha=0.3, label="", framestyle=:box)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9448600653594766,0.9448600653594766,0.9658495915032675) rtol = 0.01
        rm("concordia.png")

        # Plot single concordia line
        h = concordialine(0, 100, label="")
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.981812385620915,0.9832414705882352,0.9841176307189541) rtol = 0.01
        rm("concordia.png")
        h = concordialine(0, 100, label="", truncate=true)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.981812385620915,0.9832414705882352,0.9841176307189541) rtol = 0.01
        rm("concordia.png")

        # Plot many single concordia lines
        h = concordialine(10*randn(100).+10, 100*randn(100).+1000, label="")
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.966250588235294,0.966250588235294,0.966250588235294) rtol = 0.01
        rm("concordia.png")
        h = concordialine(10*randn(100).+10, 100*randn(100).+1000, label="", truncate=true)
        concordiacurve!(h)
        savefig(h, "concordia.png")
        img = load("concordia.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.966250588235294,0.966250588235294,0.966250588235294) rtol = 0.01
        rm("concordia.png")

        # Tera-Wasserburg analyses
        h = plot(twanalysis.(analyses), color=:blue, alpha=0.3, label="")
        twconcordiacurve!(h)
        @test h isa Plots.Plot
        savefig(h, "concordia.png")
        rm("concordia.png")

        # Tera-Wasserburg concordia
        h = twconcordialine(500, 1000, label="")
        twconcordiacurve!(h, timerange=(300, 1500))
        savefig(h, "concordia.png")
        @test load("concordia.png") isa AbstractMatrix
        rm("concordia.png")
        h = twconcordialine(500, 1000, label="", truncate=true)
        twconcordiacurve!(h)
        @test h isa Plots.Plot
        h = twconcordialine(10*randn(20).+500, 10*randn(20).+1000, label="")
        twconcordiacurve!(h)
        h = twconcordialine(10*randn(20).+500, 10*randn(20).+1000, label="", truncate=true)
        twconcordiacurve!(h)
        @test h isa Plots.Plot

        # Rank-order plot
        h = rankorder(1:10, 2*ones(10), label="")
        savefig(h, "rankorder.png")
        img = load("rankorder.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9842757843137254,0.9882103758169932,0.9906218464052285) rtol = 0.01
        rm("rankorder.png")

        h = rankorder((1:10) .± (2*ones(10)), label="")
        savefig(h, "rankorder.png")
        img = load("rankorder.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9842757843137254,0.9882103758169932,0.9906218464052285) rtol = 0.01
        rm("rankorder.png")

        # 3D analyses: 2D projection and ellipsoid
        a3 = [Analysis([t, 1+t/2, 2t], [0.2, 0.3, 0.4], [0.04 0.03 0.0; 0.03 0.09 0.05; 0.0 0.05 0.16]) for t in 0:3]
        @test plot(a3, dims=(1,3)) isa Plots.Plot
        @test plot!(plot(), a3[1]) isa Plots.Plot
        @test plot3d!(plot3d(), a3[1]) isa Plots.Plot
        h = plot3d!(plot3d(), a3, label="ellipsoids")
        savefig(h, "ell3.png")
        @test load("ell3.png") isa AbstractMatrix
        rm("ell3.png")

        # With weighted mean
        d = [10.0, 10.3, 9.8, 10.6, 10.1] .± 0.2
        m = Isoplot.wmeanband(d)
        @test m.μ ≈ value(first(wmean(d))) && m.σ ≈ 2stdev(first(wmean(d)))
        @test occursin("MSWD", m.label)
        h = rankorder(d, label="", wmean=true)
        @test h isa Plots.Plot
        @test length(h.series_list) > length(rankorder(d, label="").series_list)
        @test_throws ArgumentError rankorder([CI(1:10) for _ in 1:3], wmean=true)

        h = rankorder([CI(1:10) for _ in 1:10], label="")
        savefig(h, "CIa.png")
        img = load("CIa.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.973667565359477,0.9838965359477125,0.9901668300653594) rtol = 0.01
        rm("CIa.png")

        h = plot([CI(1:10) for _ in 1:10], label="")
        savefig(h, "CIb.png")
        img = load("CIb.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈  RGB{Float64}(0.9728901143790849,0.9831728594771241,0.9894754738562094) rtol = 0.01
        rm("CIb.png")

        h = plot([CI(1:10) for _ in 1:10], 1:10, label="")
        savefig(h, "CIc.png")
        img = load("CIc.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈  RGB{Float64}(0.9661086437908494,0.9805500816993463,0.9894015849673207) rtol = 0.01
        rm("CIc.png")

        h = plot([CI(1:10) for _ in 1:10], [CI(1:10) for _ in 1:10], label="")
        savefig(h, "CId.png")
        img = load("CId.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9829569444444444,0.9866916666666666,0.9890226307189541) rtol = 0.01
        rm("CId.png")

        # York fit
        x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
        y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
        yf = yorkfit(x, fill(0.25, 10), y, fill(0.25, 10))
        h = plot(yf, label="")
        savefig(h, "yf.png")
        img = load("yf.png")
        @test sum(img)/length(img) ≈ RGB{Float64}(0.8183448856209153, 0.9208928758169928, 0.9848070588235287) rtol = 0.01
        rm("yf.png")

        # McLean fit
        mf = mcleanfit(x, fill(0.25, 10), y, fill(0.25, 10))
        h = plot(mf, label="")
        savefig(h, "mf.png")
        img = load("mf.png")
        @test sum(img)/length(img) ≈ RGB{Float64}(0.8182011764705884, 0.9208347549019601, 0.9848014215686267) rtol = 0.01
        rm("mf.png")

        # 3D McLean fit, as 3D line and as 2D projections
        pts3 = [[t, 1 + t/2, 2t - 2] .+ 0.05 .* [sin(7t), cos(5t), sin(3t)] for t in range(0, 5, length=12)]
        mf3 = mcleanfit(pts3, [Matrix(Diagonal([0.05, 0.08, 0.1].^2)) for _ in pts3])
        h = scatter([p[1] for p in pts3], [p[2] for p in pts3], label="")
        @test plot!(h, mf3, dims=2, label="") isa Plots.Plot
        h = scatter3d([p[1] for p in pts3], [p[2] for p in pts3], [p[3] for p in pts3], label="")
        @test plot3d!(h, mf3, label="") isa Plots.Plot
        savefig(h, "mf3.png")
        @test load("mf3.png") isa AbstractMatrix
        rm("mf3.png")

        # Plot SIMS calibration line
        h = plot(calib)
        savefig(h, "calib.png")
        img = load("calib.png")
        # @test size(img) == (400,600)
        @test sum(img)/length(img) ≈ RGB{Float64}(0.9189778921568632, 0.9413858333333338, 0.9578708333333333) rtol = 0.01
        rm("calib.png")
    end

end

module MakieTest
    using Test, Statistics
    using LinearAlgebra: Diagonal
    using Measurements
    using Isoplot
    using ImageIO, FileIO, CairoMakie
    using ColorTypes

    import ..BaseTests: analyses

    const MakieExt = Base.get_extension(Isoplot, :MakieExt)

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

        # Tera-Wasserburg concordia
        f4 = Figure()
        ax4 = Axis(f4[1,1])
        l = concordialine!(ax4, 500, 1000, concordiaType=:TeraWasserburg, label="line")
        c = concordiacurve!(ax4, 300, 1500, concordiaType=:TeraWasserburg, label="curve")
        @test l isa MakieExt.ConcordiaLine
        @test c isa MakieExt.ConcordiaCurve
        axislegend(ax4)
        save("concordia.png",f4)
        @test size(load("concordia.png")) == (900, 1200)
        rm("concordia.png")
        @test_throws ArgumentError concordiacurve!(Axis(Figure()[1,1]), 0, 100, concordiaType=:Foo)

        # Many concordia lines, and a plot of many analyses at once
        f5 = Figure()
        ax5 = Axis(f5[1,1])
        concordialine!(ax5, 10*randn(20).+10, 100*randn(20).+1000)
        plot!(ax5, analyses, color=(:blue, 0.3))
        concordiacurve!(ax5) # Time range from axis limits (current axis)
        save("concordia.png",f5)
        @test size(load("concordia.png")) == (900, 1200)
        rm("concordia.png")

        # Tera-Wasserburg analyses
        f8 = Figure()
        ax8 = Axis(f8[1,1])
        plot!(ax8, twanalysis.(analyses), color=(:blue, 0.3))
        concordiacurve!(ax8, 100, 3000, concordiaType=:TeraWasserburg)
        twconcordiacurve!(ax8, 100, 3000)
        twconcordialine!(ax8, 500, 1000)
        save("concordia.png",f8)
        @test size(load("concordia.png")) == (900, 1200)
        rm("concordia.png")

        # York fit
        x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
        y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
        yf = yorkfit(x, fill(0.25, 10), y, fill(0.25, 10))
        f6 = Figure()
        ax6 = Axis(f6[1,1])
        scatter!(ax6, x, y)
        @test plot!(ax6, yf) isa MakieExt.LineFit
        @test plot!(ax6, yf, xlims=(0, 12)) isa MakieExt.LineFit

        # McLean fit
        mf = mcleanfit(x, fill(0.25, 10), y, fill(0.25, 10))
        f6 = Figure()
        ax6 = Axis(f6[1,1])
        scatter!(ax6, x, y)
        @test plot!(ax6, mf) isa MakieExt.LineFit
        @test plot!(ax6, mf, xlims=(0, 12)) isa MakieExt.LineFit

        # 3D McLean fit, as 3D line and as 2D projections
        pts3 = [[t, 1 + t/2, 2t - 2] .+ 0.05 .* [sin(7t), cos(5t), sin(3t)] for t in range(0, 5, length=12)]
        mf3 = mcleanfit(pts3, [Matrix(Diagonal([0.05, 0.08, 0.1].^2)) for _ in pts3])
        f9 = Figure()
        ax9 = Axis(f9[1,1])
        scatter!(ax9, [p[1] for p in pts3], [p[2] for p in pts3])
        @test plot!(ax9, mf3, dims=2) isa MakieExt.LineFit
        ax9₃ = Axis3(f9[1,2])
        scatter!(ax9₃, [p[1] for p in pts3], [p[2] for p in pts3], [p[3] for p in pts3])
        @test plot!(ax9₃, mf3) isa MakieExt.LineFit
        save("mf3.png", f9)
        @test size(load("mf3.png")) == (900, 1200)
        rm("mf3.png")

        # Rank-order plots
        f7 = Figure()
        ax7 = Axis(f7[1,1])
        @test rankorder!(ax7, (1:10) .± (2*ones(10))) isa MakieExt.RankOrderPlot
        @test rankorder!(ax7, 1:10, 2*ones(10)) isa MakieExt.RankOrderPlot
        @test rankorder!(ax7, [CI(1:10) for _ in 1:10]) isa MakieExt.RankOrderPlot
        # 3D analyses: ellipsoids, and 2D projections
        a3 = [Analysis([t, 1+t/2, 2t], [0.2, 0.3, 0.4], [0.04 0.03 0.0; 0.03 0.09 0.05; 0.0 0.05 0.16]) for t in 0:3]
        f3 = Figure(size=(600, 300))
        @test plot!(Axis3(f3[1,1]), a3) isa MakieExt.Ellipsoids
        @test plot!(Axis3(f3[1,2]), a3[1]) isa MakieExt.Ellipsoids
        @test plot!(Axis(f3[1,3]), Ellipse.(a3, 1, 3)) isa Makie.Poly
        save("ell3.png", f3)
        @test load("ell3.png") isa AbstractMatrix
        rm("ell3.png")
        @test rankorder!(ax7, [10.0, 10.3, 9.8, 10.6, 10.1] .± 0.2, wmean=true) isa MakieExt.RankOrderPlot
    end
end

