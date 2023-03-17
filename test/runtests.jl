using Isoplot
using Test, Statistics

@testset "UPb" begin
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

    @test rand(d1) isa Vector{Float64}
    @test rand(d1, 10) isa Matrix{Float64}

    x = [22.70307499779583, 22.681635852743682, 22.63876085494785, 22.61732500220417, 22.638764147256317, 22.68163914505215, 22.70307499779583]
    y = [0.4089925091091875, 0.40901969166358015, 0.4086701825543926, 0.40829349089081246, 0.4082663083364198, 0.4086158174456074, 0.4089925091091875]
    e1 = ellipse(d1, npoints=7)
    @test e1 isa Isoplot.Shape
    @test e1.x ≈ x
    @test e1.y ≈ y

    tₗₗ = 35
    ui = upper_intercept(tₗₗ, d1)
    @test ui isa Isoplot.Measurement

    N = 10000
    uis = upper_intercept(tₗₗ, d1, N)
    @test uis isa Vector{Float64}
    @test mean(uis) ≈ ui.val atol=(4*ui.err/sqrt(N))
    @test std(uis) ≈ ui.err rtol=0.03

end
