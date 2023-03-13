
# Make an ellipse from a UPbAnalysis object
function ellipse(d::UPbAnalysis;
        sigmalevel=sqrt(invlogccdf(Chisq(2), log(0.05))),
        npoints::Integer=100,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return ellipse(d.μ[1], d.μ[2], a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function ellipse(x₀, y₀, a, b, θ; npoints::Integer=100)
    t = range(0, 2π, length=npoints)
    x = a*cos(θ)*cos.(t) .- b*sin(θ)*sin.(t) .+ x₀
    y = a*sin(θ)*cos.(t) .+ b*cos(θ)*sin.(t) .+ y₀
    return Shape(x, y)
end
export ellipse

# Non-exported function: return semimajor and minor axes for a given U-Pb analysis
function ellipseparameters(d::UPbAnalysis{T}, sigmalevel::Number) where T
    if !any(isnan, d.Σ)
        # Calculate eigenvectors and eigenvalues from the covariance matrix.
        # V: matrix of eigenvectors, D: diagonal matrix of eigenvalues
        F = eigen(d.Σ)
        # Find index of semimajor and semiminor axes
        major = argmax(F.values)
        minor = argmin(F.values)
        v = view(F.vectors, :, major)

        # Calculate angle of major axis of ellipse from horizontal
        θ = atan(v[2]/v[1])

        # Calculate length of semimajor and semiminor axes for given p-value
        a = T(sigmalevel)*sqrt(abs(F.values[major]))
        b = T(sigmalevel)*sqrt(abs(F.values[minor]))

        return a, b, θ
    else
        return T.((NaN, NaN, NaN))
    end
end
