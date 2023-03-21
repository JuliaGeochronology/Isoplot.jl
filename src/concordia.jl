
Δ68(t,(slope,r75,r68)) = slope * (exp(λ235U.val*t) - 1 - r75) + r68 - exp(λ238U.val*t) + 1
dΔ68(t,(slope,r75,r68)) = slope * λ235U.val*exp(λ235U.val*t) - λ238U.val*exp(λ238U.val*t)

function newton_zero(f, df, x0, args::Tuple, iterations=10)
    for i in 1:iterations
        δx = f(x0, args)/df(x0, args)
        x0 -= δx
    end
    return x0
end

function upperintercept(tₗₗ::Number, s::Shape{T,T};
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
    ) where T

    # Get ratios from our ellipse
    r75, r68 = s.x, s.y
    r75₀, r68₀ = center(s)
    # Return early if our lead loss time is too old or anything is NaN'd
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
    0 < slope₋ < Inf || return T(NaN) ± T(NaN)
    r75₊, r68₊ = argmin(x->atan((x[2]-r68ₗₗ)/(x[1]-r75ₗₗ)), zip(r75, r68))
    slope₊ = (r68₊-r68ₗₗ)/(r75₊-r75ₗₗ)
    0 < slope₊ < Inf || return T(NaN) ± T(NaN)

    ui₀ = newton_zero(Δ68, dΔ68, 4.567e3, (slope₀,r75₀,r68₀))
    # Return early if our upper intercept is younger than the analysis
    ui₀ > log(r68₀+1)/λ238U.val || return T(NaN) ± T(NaN)
    ui₋ = newton_zero(Δ68, dΔ68, 4.567e3, (slope₋,r75₋,r68₋))
    ui₊ = newton_zero(Δ68, dΔ68, 4.567e3, (slope₊,r75₊,r68₊))

    return ui₀ ± (ui₊ - ui₋)/2sigmalevel
end

function upperintercept(tₗₗ::Number, d::UPbAnalysis{T}, nresamplings::Integer) where T
    ui = zeros(T, nresamplings)

    # Get ratios
    r75₀, r68₀ = d.μ
    # Return early if our lead loss time is too old or anything is NaN'd
    tₗₗ < log(r68₀+1)/λ238U.val || return fill!(uis, T(NaN))
    tₗₗ < log(r75₀+1)/λ235U.val || return fill!(uis, T(NaN))

    # Calculate isotopic ratios of our time of Pb-loss
    r75ₗₗ = exp(λ235U.val*tₗₗ) - 1
    r68ₗₗ = exp(λ238U.val*tₗₗ) - 1
    slope₀ = (r68₀-r68ₗₗ)/(r75₀-r75ₗₗ)

    samples = rand(d, nresamplings)
    for i in axes(samples,2)
        r75, r68 = view(samples, :, i)
        slope = (r68-r68ₗₗ)/(r75-r75ₗₗ)
        ui[i] = newton_zero(Δ68, dΔ68, 4.567e3, (slope,r75,r68))
    end
    return ui
end

upperintercept(tₗₗ::Number, d::UPbAnalysis) = upperintercept(tₗₗ, ellipse(d; npoints=50))

function upperintercept(d::Vector{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    ui = zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    for i in eachindex(slopes, intercepts)
        ui[i] = newton_zero(Δ68, dΔ68, 4.567e3, (slopes[i],zero(T),intercepts[i]))
    end
    return ui
end

function lowerintercept(d::Vector{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    li = zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    for i in eachindex(slopes, intercepts)
        li[i] = newton_zero(Δ68, dΔ68, 0.0, (slopes[i],zero(T),intercepts[i]))
    end
    return li
end

function intercepts(d::Vector{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    ui, li = zeros(T, nresamplings), zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    for i in eachindex(slopes, intercepts)
        ui[i] = newton_zero(Δ68, dΔ68, 4.567e3, (slopes[i],zero(T),intercepts[i]))
        li[i] = newton_zero(Δ68, dΔ68, 0.0, (slopes[i],zero(T),intercepts[i]))
    end
    return ui, li
end

function fit_lines(d::Vector{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    nanalyses = length(d)
    # Vector of ratios
    r68 = zeros(T, nanalyses)
    # Draw random ratios from each analysis
    randratios = rand.(d, nresamplings)
    # Allocate temporary arrays for regression
    A = ones(T, nanalyses, 2)
    # Allocate output slopes and intercepts
    slopes, intercepts = zeros(T, nresamplings), zeros(T, nresamplings)
    @inbounds for n in eachindex(slopes, intercepts)
        for i in eachindex(d, randratios)
            A[i,2] = randratios[i][1,n]
            r68[i] = randratios[i][2,n]
        end
        # Linear regression
        ϕ = A\r68 # or perhgaps alternatively in some cases ϕ = ldiv!(lu!(A), r68)?
        slopes[n], intercepts[n] = ϕ[2], ϕ[1]
    end
    return slopes, intercepts
end
