
Δ68(t,(slope,r75,r68)) = slope * (exp(value(λ235U)*t) - 1 - r75) + r68 - exp(value(λ238U)*t) + 1
dΔ68(t,(slope,r75,r68)) = slope * value(λ235U)*exp(value(λ235U)*t) - value(λ238U)*exp(value(λ238U)*t)

function newton_zero(f, df, x0, args::Tuple, iterations=10)
    for i in 1:iterations
        δx = f(x0, args)/df(x0, args)
        x0 -= δx
    end
    return x0
end

# Pb207/Pb206 age
age76(d::UPbAnalysis{T}) where {T} = upperintercept(zero(T), d)
age76(d::UThPbAnalysis) = age76(UPbAnalysis(d))

# 2D Concordia age (via probability density of analysis intersecting the 2d concordia line)
function ageconcordia(d::UPbAnalysis)
    isnan(d) && return NaN ± NaN
    isposdef(cov(d)) || return NaN ± NaN
    dist = MvNormal(mean(d), cov(d))
    μ, mswd = wmean(SVector(age75(d), age68(d)), corrected=true)

    # Evaluate at 256 points between -6σ and  +6σ
    agerange = range(value(μ)-6stdev(μ), value(μ)+6stdev(μ), length=256)
    intersection = @SVector [pdf(dist, SVector(ratio(agerange[i], value(λ235U)), ratio(agerange[i], value(λ238U)))) for i in 1:256]
    
    # Return Measurement
    return histmean(intersection, agerange) ± histstd(intersection, agerange, corrected=false)
end
ageconcordia(d::UThPbAnalysis) = ageconcordia(UPbAnalysis(d))

# 3-D Concordia age (via probability density of analysis intersecting the 3d concordia line)
function ageconcordia3d(d::UThPbAnalysis)
    isnan(d) && return NaN ± NaN
    isposdef(cov(d)) || return NaN ± NaN
    dist = MvNormal(mean(d), cov(d))
    μ, mswd = wmean(SVector(age75(d), age68(d), age82(d)), corrected=true)

    # Evaluate at 256 points between -6σ and  +6σ
    agerange = range(value(μ)-6stdev(μ), value(μ)+6stdev(μ), length=256)
    intersection = @SVector [pdf(dist, SVector(ratio(agerange[i], value(λ235U)), ratio(agerange[i], value(λ238U)), ratio(agerange[i], value(λ232Th)))) for i in 1:256]
    
    # Return Measurement
    return histmean(intersection, agerange) ± histstd(intersection, agerange, corrected=false)
end

function isconcordant(d::UPbAnalysis;
        sigmalevel::Number=2.447746830680816, # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))
        npoints::Integer = 50,
    )
    isnan(d) && return false

    # Make error ellipse
    e = Ellipse(d; sigmalevel, npoints)
    # Calculate mean
    μ, mswd = wmean([age75(d), age68(d)], corrected=true)

    # Evaluate at 256 points between -6σ and +6σ
    agerange = range(value(μ)-6stdev(μ), value(μ)+6stdev(μ), length=256)
    for i in eachindex(agerange)
        xᵢ = ratio(agerange[i], value(λ235U))
        yᵢ = ratio(agerange[i], value(λ238U))
        if inpolygon(e.x, e.y, (xᵢ, yᵢ))
            return true
        end
    end
    return false
end
isconcordant(d::UThPbAnalysis) = isconcordant(UPbAnalysis(d))

function upperintercept(tₗₗ::Number, s::Ellipse{T}, sigmalevel::T=2.447746830680816) where {T<:AbstractFloat}
    # bivariate p=0.05 level: sqrt(invlogccdf(Chisq(2), log(0.05)))

    # Get ratios from our ellipse
    r75, r68 = s.x, s.y
    r75₀, r68₀ = s.x₀, s.y₀
    σ75₀, σ68₀ = s.σx₀, s.σy₀
    # Return early if our lead loss time is too old or anything is NaN'd
    tₗₗ < age(r68₀,value(λ238U)) || return T(NaN) ± T(NaN)
    tₗₗ < age(r75₀,value(λ235U)) || return T(NaN) ± T(NaN)

    # If reversely discordant, move to the closest point on Concordia rather
    # than projecting down a fictive "lead gain" array, increasing uncertainty
    # by sqrt(MSWD) if discordance is large
    age68 = age(r68₀ ± σ68₀, value(λ238U))
    age75 = age(r75₀ ± σ75₀, value(λ235U))
    age75.val > age68.val || return first(wmean([age68, age75], corrected=true))

    # Calculate isotopic ratios of our time of Pb-loss
    r75ₗₗ = ratio(tₗₗ, value(λ235U))
    r68ₗₗ = ratio(tₗₗ, value(λ238U))
    slope₀ = (r68₀-r68ₗₗ)/(r75₀-r75ₗₗ)

    # Find the values on the margin of the ellipse with the
    # largest and smallest angular difference from the center
    r75₋, r68₋ = argmax(x->atan((x[2]-r68ₗₗ)/(x[1]-r75ₗₗ)), zip(r75, r68))
    slope₋ = (r68₋-r68ₗₗ)/(r75₋-r75ₗₗ)
    0 < slope₋ < Inf || return T(NaN) ± T(NaN)
    r75₊, r68₊ = argmin(x->atan((x[2]-r68ₗₗ)/(x[1]-r75ₗₗ)), zip(r75, r68))
    slope₊ = (r68₊-r68ₗₗ)/(r75₊-r75ₗₗ)
    0 < slope₊ < Inf || return T(NaN) ± T(NaN)

    # Find the upper intercept of our Pb-loss arrays with Concordia
    ui₀ = newton_zero(Δ68, dΔ68, t🜨, (slope₀,r75₀,r68₀))
    0 < ui₀ < t🜨 || return T(NaN) ± T(NaN)
    ui₋ = newton_zero(Δ68, dΔ68, t🜨, (slope₋,r75₋,r68₋))
    ui₊ = newton_zero(Δ68, dΔ68, t🜨, (slope₊,r75₊,r68₊))
    # Direct uncertainty, from spread in intercepts given size of ellipse
    σ = (value(ui₊) - value(ui₋))/2sigmalevel
    # Include also uncertainty, from lower intercept if tₗₗ (and ui) are `Measurement`s
    return value(ui₀) ± σcombined(ui₀, σ)
end
σcombined(m::Measurement, σ) = sqrt(stdev(m)^2 + σ^2)
σcombined(m, σ) = σ # If m is not a Measurement

upperintercept(tₗₗ::Number, d::UPbAnalysis) = upperintercept(tₗₗ, Ellipse(d; npoints=50))

function upperintercept(tₗₗ::Number, d::UPbAnalysis{T}, nresamplings::Integer) where T
    # Get ratios
    r75₀, r68₀ = mean(d)
    # Return early if our lead loss time is too old or anything is NaN'd
    tₗₗ < log(r68₀+1)/value(λ238U) || return fill(T(NaN), nresamplings)
    tₗₗ < log(r75₀+1)/value(λ235U) || return fill(T(NaN), nresamplings)
    
    # Calculate isotopic ratios of our time of Pb-loss
    r75ₗₗ = exp(value(λ235U)*tₗₗ) - 1
    r68ₗₗ = exp(value(λ238U)*tₗₗ) - 1
    slope₀ = (r68₀-r68ₗₗ)/(r75₀-r75ₗₗ)

    ui = zeros(T, nresamplings)
    samples = rand(d, nresamplings)
    @assert axes(samples,2) == eachindex(ui)
    @inbounds for i in axes(samples,2)
        r75, r68 = view(samples, :, i)
        slope = (r68-r68ₗₗ)/(r75-r75ₗₗ)
        ui[i] = newton_zero(Δ68, dΔ68, t🜨, (slope,r75,r68))
    end
    return ui
end

function upperintercept(d::Collection{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    ui = zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    @inbounds for i in eachindex(ui, slopes, intercepts)
        ui[i] = newton_zero(Δ68, dΔ68, t🜨, (slopes[i],zero(T),intercepts[i]))
    end
    return ui
end
function upperintercept(d::Collection{UPbAnalysis{T}}) where {T}
    yf = yorkfit(d)
    return newton_zero(Δ68, dΔ68, t🜨, (yf.slope, yf.xm, yf.ym))
end

function lowerintercept(d::Collection{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    li = zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    @inbounds for i in eachindex(li, slopes, intercepts)
        li[i] = newton_zero(Δ68, dΔ68, zero(T), (slopes[i],zero(T),intercepts[i]))
    end
    return li
end
function lowerintercept(d::Collection{UPbAnalysis{T}}) where {T}
    yf = yorkfit(d)
    return newton_zero(Δ68, dΔ68, zero(T), (yf.slope, yf.xm, yf.ym))
end

function intercepts(d::Collection{UPbAnalysis{T}}, nresamplings::Integer) where {T}
    ui, li = zeros(T, nresamplings), zeros(T, nresamplings)
    slopes, intercepts = fit_lines(d, nresamplings)
    @inbounds for i in eachindex(ui, li, slopes, intercepts)
        ui[i] = newton_zero(Δ68, dΔ68, t🜨, (slopes[i],zero(T),intercepts[i]))
        li[i] = newton_zero(Δ68, dΔ68, zero(T), (slopes[i],zero(T),intercepts[i]))
    end
    return ui, li
end
function intercepts(d::Collection{UPbAnalysis{T}}) where {T}
    yf = yorkfit(d)
    ui = newton_zero(Δ68, dΔ68, t🜨, (yf.slope, yf.xm, yf.ym))
    li = newton_zero(Δ68, dΔ68, zero(T), (yf.slope, yf.xm, yf.ym))
    return ui, li
end

function fit_lines(d::Collection{UPbAnalysis{T}}, nresamplings::Integer) where {T}
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
