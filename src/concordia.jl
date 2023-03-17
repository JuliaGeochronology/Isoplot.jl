
# Make an ellipse from a UPbAnalysis object
function ellipse(d::UPbAnalysis;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer=50,
    )
    a, b, θ = ellipseparameters(d, sigmalevel)
    return ellipse(d.μ[1], d.μ[2], a, b, θ; npoints)
end
# Make an ellipse if given x and y positions, major and minor axes, and rotation
function ellipse(x₀, y₀, a, b, θ; npoints::Integer=50)
    t = range(0, 2π, length=npoints)
    x = a*cos(θ)*cos.(t) .- b*sin(θ)*sin.(t) .+ x₀
    y = a*sin(θ)*cos.(t) .+ b*cos(θ)*sin.(t) .+ y₀
    return Shape(x, y)
end
export ellipse

# Non-exported function: return semimajor and minor axes for a given U-Pb analysis
function ellipseparameters(d::UPbAnalysis{T}, sigmalevel::Number) where T

    # Quickly exit if any NaNs
    any(isnan, d.Σ) && return T.((NaN, NaN, NaN))

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

end

Δ68(t,slope,r75,r68) = slope * (exp(λ235U.val*t) - 1 - r75) + r68 - exp(λ238U.val*t) + 1

function upper_intercept(tₗₗ::Number, s::Shape{T,T};
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
    ) where T

    # Get ratios from our ellipse
    r75, r68 = s.x, s.y
    r75₀, r68₀ = center(s)
    # Return early if our lead loss time is too young or anything is NaN'd
    tₗₗ < log(r68₀+1)/λ238U.val || return T(NaN) ± T(NaN)
    tₗₗ < log(r75₀+1)/λ235U.val || return T(NaN) ± T(NaN)

    # Calculate isotopic ratios of our time of Pb-loss
    r75ₗₗ = exp(λ235U.val*tₗₗ) - 1
    r68ₗₗ = exp(λ238U.val*tₗₗ) - 1
    slope₀ = (r68₀-r68ₗₗ)/(r75₀-r75ₗₗ)

    # Find the values on the margin of the ellipse with the
    # largest and smallest angular difference from the center
    r75₋, r68₋ = argmax(x->atan((x[2]-r68ₗₗ)/(x[1]-r75ₗₗ)), zip(r75, r68))
    slope₋ = (r68₋-r68ₗₗ)/(r75₋-r75ₗₗ)
    r75₊, r68₊ = argmin(x->atan((x[2]-r68ₗₗ)/(x[1]-r75ₗₗ)), zip(r75, r68))
    slope₊ = (r68₊-r68ₗₗ)/(r75₊-r75ₗₗ)

    ui₀ = find_zero(t->Δ68(t,slope₀,r75₀,r68₀), 4.567e3)
    ui₋ = find_zero(t->Δ68(t,slope₋,r75₋,r68₋), 4.567e3)
    ui₊ = find_zero(t->Δ68(t,slope₊,r75₊,r68₊), 4.567e3)

    return ui₀ ± (ui₊ - ui₋)/2sigmalevel
end

function upper_intercept(tₗₗ::Number, d::UPbAnalysis{T}, nresamplings::Integer) where T
    results = zeros(T, nresamplings)

    # Get ratios
    r75₀, r68₀ = d.μ
    # Return early if our lead loss time is too young or anything is NaN'd
    tₗₗ < log(r68₀+1)/λ238U.val || return fill!(results, T(NaN))
    tₗₗ < log(r75₀+1)/λ235U.val || return fill!(results, T(NaN))

    # Calculate isotopic ratios of our time of Pb-loss
    r75ₗₗ = exp(λ235U.val*tₗₗ) - 1
    r68ₗₗ = exp(λ238U.val*tₗₗ) - 1
    slope₀ = (r68₀-r68ₗₗ)/(r75₀-r75ₗₗ)

    samples = rand(d, nresamplings)
    for i in axes(samples,2)
        r75, r68 = view(samples, :, i)
        slope = (r68-r68ₗₗ)/(r75-r75ₗₗ)
        ui = find_zero(t->Δ68(t,slope,r75,r68), 4.567e3)
        results[i] = ui
    end
    return results
end

upper_intercept(tₗₗ::Number, d::UPbAnalysis) = upper_intercept(tₗₗ, ellipse(d; npoints=100))

export upper_intercept
