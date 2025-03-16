
"""
```julia
dist_ll(dist::Collection, mu::Collection, sigma::Collection, tmin::Number, tmax::Number)
dist_ll(dist::Collection, analyses::Collection{<:Measurement}, tmin::Number, tmax::Number)
```
Return the log-likelihood of a set of mineral ages with means `mu` and
uncertianty `sigma` being drawn from a given source (i.e., crystallization / closure)
distribution `dist`, with terms to prevent runaway at low N.

### Examples
```julia
mu, sigma = collect(100:0.1:101), 0.01*ones(11)
ll = dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100, 101)
```
"""
function dist_ll(dist::Collection, mu::Collection{<:Number}, sigma::Collection{<:Number}, tmin::Number, tmax::Number)
    tmax >= tmin || return NaN
    any(isnan, mu) && return NaN
    any(x->!(x>0), sigma) && return NaN
    nbins = length(dist) - 1
    Δt = abs(tmax-tmin)
    dt = Δt/nbins

    # Cycle through each datum in dataset
    ll = zero(float(eltype(dist)))
    @inbounds for j in eachindex(mu, sigma)
        μⱼ, σⱼ = mu[j], sigma[j]

        # Find equivalent index position of μⱼ in the `dist` array
        ix = (μⱼ - tmin) / Δt * nbins + 1
        # If possible, prevent aliasing problems by interpolation
        if (σⱼ < Δt / nbins) && 1 < ix < length(dist)
            # Interpolate corresponding distribution value
            f = floor(Int, ix)
            δ = ix - f
            likelihood = (dist[f+1]*δ + dist[f]*(1-δ)) / Δt
        else
            # Otherwise, sum contributions from Gaussians at each point in distribution
            𝑖 = 1:length(dist)
            likelihood = zero(float(eltype(dist)))
            normconst = 1/(length(dist) * σⱼ * sqrt(2 * pi))
            @turbo for i in eachindex(dist, 𝑖)
                distx = tmin + dt * (𝑖[i] - 1) # time-position of distribution point
                # Likelihood curve follows a Gaussian PDF. Note: Δt cancels
                likelihood += dist[i] * normconst * exp(-(distx - μⱼ)^2 / (2 * σⱼ * σⱼ))
            end
        end
        ll += log(likelihood)
    end
    return ll
end
function dist_ll(dist::Collection, analyses::Collection{<:Measurement}, tmin::Number, tmax::Number)
    tmax >= tmin || return NaN
    any(isnan, analyses) && return NaN
    any(x->!(stdev(x) > 0), analyses) && return NaN
    nbins = length(dist) - 1
    Δt = abs(tmax - tmin)
    dt = Δt/nbins

    # Cycle through each datum in dataset
    ll = zero(float(eltype(dist)))
    @inbounds for j in eachindex(analyses)
        dⱼ = analyses[j]
        μⱼ, σⱼ = value(dⱼ), stdev(dⱼ)

        # Find equivalent index position of μⱼ in the `dist` array
        ix = (μⱼ - tmin) / Δt * nbins + 1
        # If possible, prevent aliasing problems by interpolation
        if (σⱼ < Δt / nbins) && 1 < ix < length(dist)
            # Interpolate corresponding distribution value
            f = floor(Int, ix)
            δ = ix - f
            likelihood = (dist[f+1]*δ + dist[f]*(1-δ)) / Δt
        else
            # Otherwise, sum contributions from Gaussians at each point in distribution
            𝑖 = 1:length(dist)
            likelihood = zero(float(eltype(dist)))
            normconst = 1/(length(dist) * σⱼ * sqrt(2 * pi))
            @inbounds @fastmath @simd ivdep for i in eachindex(dist, 𝑖)
                distx = tmin + dt * (𝑖[i] - 1) # time-position of distribution point
                # Likelihood curve follows a Gaussian PDF. Note: Δt cancels
                likelihood += dist[i] * normconst * exp(-(distx - μⱼ)^2 / (2 * σⱼ * σⱼ))
            end
        end
        ll += log(likelihood)
    end
    return ll
end
function dist_ll(dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}, tmin::Number, tmax::Number, tll::Number) where {T<:AbstractFloat}
    tmax >= tmin || return T(NaN)
    any(isnan, analyses) && return T(NaN)

    tbinedges = range(tmin, tmax, length=length(dist))
    @assert eachindex(tbinedges) == eachindex(dist)
    dt = step(tbinedges)
    Σdist = sum(dist)
    r75ₗₗ = ratio(tll, value(λ235U))
    r68ₗₗ = ratio(tll, value(λ238U))
    μₗₗ = SVector(r75ₗₗ, r68ₗₗ)

    # Cycle through each datum in dataset
    ll = zero(float(eltype(dist)))
    @inbounds for j in eachindex(analyses)
        d = analyses[j]
        μⱼ = mean(d)
        Σⱼ = cov(d)

        # Cycle through each time step
        likelihood = zero(T)
        r75ᵢ = ratio(first(tbinedges)-step(tbinedges), value(λ235U))
        r68ᵢ = ratio(first(tbinedges)-step(tbinedges), value(λ238U))
        for i in eachindex(dist)
            μlast = SVector(r75ᵢ, r68ᵢ)

            # Rotation matrix that would rotate discordia line between tll and tbinedges[i] to vertical
            r75ᵢ = ratio(tbinedges[i], value(λ235U))
            r68ᵢ = ratio(tbinedges[i], value(λ238U))
            R = RotMatrix(π/2 - atan(r68ᵢ-r68ₗₗ, r75ᵢ-r75ₗₗ))

            # Rotate means and covariance matrix, with proposed time of Pb-loss at origin
            μlastᵣ = R * (μlast-μₗₗ)
            μᵣ = R * (μⱼ-μₗₗ)
            Σᵣ = R * Σⱼ * R'
            μ₁, σ₁ = first(μᵣ), sqrt(first(Σᵣ))

            # Product of PDF of marginal disribution of rotated bivariate Gaussian and `dist`
            dμ₁ = abs(first(μlastᵣ)*last(μᵣ)/last(μlastᵣ))
            likelihood += pdf(Normal(μ₁, σ₁), zero(T))*dμ₁*dist[i]/(dt*Σdist)
        end
        ll += log(likelihood)
    end
    return ll
end


# Encoding of additional prior assumptions about likely age given observed dispersion
function prior_ll(mu::Collection, sigma::Collection, tmin::Number, tmax::Number)
    i₋ = argmin(mu)
    mu₋, sigma₋ = mu[i₋], sigma[i₋]
    i₊ = argmax(mu)
    mu₊, sigma₊ = mu[i₊], sigma[i₊]
    # Calculate a weighted mean and examine our MSWD
    (wm, wsigma, mswd) = wmean(mu, sigma, corrected=true) 
    # Degrees of freedom
    f = length(mu) - 1
    return prior_ll(wm, wsigma, mswd, mu₋, sigma₋, mu₊, sigma₊, f, tmin, tmax)
end
function prior_ll(analyses::Collection{<:Measurement}, tmin::Number, tmax::Number)
    a₋ = minimum(analyses)
    a₊ = maximum(analyses)
    # Calculate a weighted mean and examine our MSWD
    (wm, mswd) = wmean(analyses, corrected=true)
    # Degrees of freedom
    f = length(analyses) - 1
    return prior_ll(wm.val, wm.err, mswd, a₋.val, a₋.err, a₊.val, a₊.err, f, tmin, tmax)
end
function prior_ll(wm, wsigma, mswd, mu₋, sigma₋, mu₊, sigma₊, f, tmin, tmax)
    # Height of MSWD distribution relative to height at MSWD = 1
    # (see Wendt and Carl, 1991, Chemical Geology)
    Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1)) * (f > 0)
    Zf = max(min(Zf, 1.0), 0.0)
    # To avoid instability of the MCMC for small datasets (low N),
    # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
    # the youngest-zircon interpretation at low Zf (MSWD far from one). The
    # specific penalty factors applied here were determined by training  
    # against synthetic datasets in Keller et al. 2018.
    # In other words, these are basically context-dependent prior distributions on tmax and tmin
    ll = -(2/log(2+f)) * (                                # Scaling factor that decreases with log number of data points (i.e., no penalty at high N)
        log((abs(tmin - wm)+wsigma)/wsigma)*Zf +          # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmax - wm)+wsigma)/wsigma)*Zf +          # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmin - mu₋)+sigma₋)/sigma₋)*(1-Zf) +     # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
        log((abs(tmax - mu₊)+sigma₊)/sigma₊)*(1-Zf)       # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    )
    return ll
end


"""
```julia
metropolis_min(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; burnin::Integer=0, t0prior=Uniform(0,minimum(value.(age68.(analyses)))), lossprior=Uniform(0,100))
metropolis_min(nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)
metropolis_min(nsteps::Integer, dist::Collection, analyses::Collection{<:UPbAnalysis; burnin::Integer=0)
```
Run a Metropolis sampler to estimate the minimum of a finite-range source
distribution `dist` using samples drawn from that distribution -- e.g., estimate
zircon eruption ages from a distribution of zircon crystallization ages.

### Examples
```julia
tmindist = metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)

tmindist, t0dist = metropolis_min(2*10^5, HalfNormalDistribution, analyses, burnin=10^5)
```
"""
metropolis_min(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; kwargs...) = metropolis_min(nsteps, dist, value.(data), stdev.(data); kwargs...)
function metropolis_min(nsteps::Integer, dist::Collection, mu::Collection, sigma::Collection; kwargs...)
    # Allocate ouput array
    tmindist = Array{float(eltype(mu))}(undef,nsteps)
    # Run Metropolis sampler
    return metropolis_min!(tmindist, dist, mu, sigma; kwargs...)
end
function metropolis_min(nsteps::Integer, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; kwargs...) where {T}
    # Allocate ouput arrays
    tmindist = Array{T}(undef,nsteps)
    t0dist = Array{T}(undef,nsteps)
    # Run Metropolis sampler
    metropolis_min!(tmindist, t0dist, dist, analyses; kwargs...)
    return tmindist, t0dist
end


"""
```julia
metropolis_min!(tmindist::DenseArray, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)
metropolis_min!(tmindist::DenseArray, t0dist::DenseArray, dist::Collection, analyses::Collection{<:UPbAnalysis}; burnin::Integer=0) where {T}
```
In-place (non-allocating) version of `metropolis_min`, fills existing array `tmindist`.

Run a Metropolis sampler to estimate the minimum of a finite-range source
distribution `dist` using samples drawn from that distribution -- e.g., estimate
zircon eruption ages from a distribution of zircon crystallization ages.

### Examples
```julia
metropolis_min!(tmindist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
```
"""
function metropolis_min!(tmindist::DenseArray{<:Number}, dist::Collection{<:Number}, mu::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}; burnin::Integer=0)

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9

    # Initial step sigma for Gaussian proposal distributions
    youngest, oldest = minimum(mu), maximum(mu)
    Δt = (oldest - youngest) + sqrt(sigma[argmin(mu)]^2 + sigma[argmax(mu)]^2)
    tminstep = tmaxstep = Δt / length(mu)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = youngest - sigma[argmin(mu)]
    tmaxₚ = tmax = oldest + sigma[argmax(mu)]

    # Log likelihood of initial proposal
    llₚ = ll = dist_ll(dist, mu, sigma, tmin, tmax) + prior_ll(mu, sigma, tmin, tmax)

    # Burnin
    for i=1:burnin
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tminstep*randn())
        (r > 0.5) && (tminₚ += tmaxstep*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu, sigma, tminₚ, tmaxₚ) + prior_ll(mu, sigma, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ-tmax)*stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist)
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tminstep*randn())
        (r > 0.5) && (tminₚ += tmaxstep*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu, sigma, tminₚ, tmaxₚ) + prior_ll(mu, sigma, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ-tmax)*stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
        end
        tmindist[i] = tmin
    end
    return tmindist
end
function metropolis_min!(tmindist::DenseArray{T}, tlldist::DenseArray{T}, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; burnin::Integer=0, tllprior=Uniform(0,minimum(value.(age68.(analyses)))), lossprior=Uniform(0,100), method=:projection) where {T<:AbstractFloat}
    @assert eachindex(tmindist) == eachindex(tlldist)
    @assert (method === :projection || method === :bivariate) "Allowed methods are `:projection` or `:bivariate`"

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9

    # Process input analyses
    ellipses = Ellipse.(analyses)
    ages68 = log.(one(T) .+ (ellipses .|> e->e.y₀))./value(λ238U)
    ages = @. upperintercept(zero(T), ellipses)

    # Initial step sigma for Gaussian proposal distributions
    youngest, oldest = minimum(ages), maximum(ages)
    Δt = (oldest.val - youngest.val) + sqrt(+ oldest.err^2 + youngest.err^2)
    tminstep = tmaxstep = Δt / length(analyses)
    tllstep = youngest.val/50

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = value(youngest)
    tmaxₚ = tmax = value(oldest)
    tllₚ = tll = zero(T)

    # Ensure priors for time and amount of lead loss are appropriately truncated
    tllprior = truncated(tllprior, 0, minimum(ages68))
    lossprior = truncated(lossprior, 0, 100) # percent

    # Log likelihood of initial proposal
    ll = logpdf(tllprior, tll)
    if method === :projection
        ll += dist_ll(dist, ages, tmin, tmax) 
    elseif method === :bivariate
        ll += dist_ll(dist, analyses, tmin, tmax, tll) 
    end
    ll += prior_ll(ages, tmin, tmax)
    for i in eachindex(ages, ages68)
        loss = 100*max(one(T) - (ages68[i] - tll) / (value(ages[i]) - tll), zero(T))
        ll += logpdf(lossprior, loss)
    end
    llₚ = ll

    # Burnin
    for i = 1:burnin
        tminₚ, tmaxₚ, tllₚ = tmin, tmax, tll
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tminstep * randn()
        elseif r < 0.70
            tmaxₚ += tmaxstep * randn()
        else
            tllₚ += tllstep * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(tllₚ, ellipses)
        llₚ = logpdf(tllprior, tllₚ)
        if method === :projection
            llₚ += dist_ll(dist, ages, tminₚ, tmaxₚ)
        elseif method === :bivariate
            llₚ += dist_ll(dist, analyses, tminₚ, tmaxₚ, tllₚ)
        end
        llₚ += prior_ll(ages, tminₚ, tmaxₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - tllₚ) / (value(ages[i]) - tllₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            tll = tllₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist, tlldist)
        tminₚ, tmaxₚ, tllₚ = tmin, tmax, tll
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tminstep * randn()
        elseif r < 0.70
            tmaxₚ += tmaxstep * randn()
        else
            tllₚ += tllstep * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(tllₚ, ellipses)
        llₚ = logpdf(tllprior, tllₚ)
        if method === :projection
            llₚ += dist_ll(dist, ages, tminₚ, tmaxₚ)
        elseif method === :bivariate
            llₚ += dist_ll(dist, analyses, tminₚ, tmaxₚ, tllₚ)
        end
        llₚ += prior_ll(ages, tminₚ, tmaxₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - tllₚ) / (value(ages[i]) - tllₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            tll = tllₚ
        end
        tmindist[i] = tmin
        tlldist[i] = tll
    end
    return tmindist
end

"""
```julia
metropolis_minmax(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; burnin::Integer=0)
metropolis_minmax(nsteps::Integer, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)
```
Run a Metropolis sampler to estimate the extrema of a finite-range source
distribution `dist` using samples drawn from that distribution -- e.g.,
estimate zircon saturation and eruption ages from a distribution of zircon
crystallization ages.

### Examples
```julia
tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
```
"""
metropolis_minmax(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; kwargs...) = metropolis_minmax(nsteps, dist, value.(data), stdev.(data); kwargs...)
function metropolis_minmax(nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; kwargs...)
    # Allocate ouput arrays
    acceptancedist = falses(nsteps)
    lldist = Array{float(eltype(dist))}(undef,nsteps)
    tmaxdist = Array{float(eltype(mu))}(undef,nsteps)
    tmindist = Array{float(eltype(mu))}(undef,nsteps)
    # Run metropolis sampler
    return metropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, dist, mu, sigma; kwargs...)
end
function metropolis_minmax(nsteps::Integer, dist::Collection{T}, analyses::Collection{<:UPbAnalysis{T}}; kwargs...) where {T<:AbstractFloat}
    # Allocate ouput arrays
    acceptancedist = falses(nsteps)
    lldist = Array{T}(undef,nsteps)
    t0dist = Array{T}(undef,nsteps)
    tmaxdist = Array{T}(undef,nsteps)
    tmindist = Array{T}(undef,nsteps)
    # Run metropolis sampler
    return metropolis_minmax!(tmindist, tmaxdist, t0dist, lldist, acceptancedist, dist, analyses; kwargs...)
end

"""
```julia
metropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, nsteps::Integer, dist::AbstractArray, data::AbstractArray, uncert::AbstractArray; burnin::Integer=0)
metropolis_minmax!(tmindist, tmaxdist, t0dist, lldist, acceptancedist, nsteps::Integer, dist::Collection, analyses::Collection{<:UPbAnalysis}; burnin::Integer=0)
```
In-place (non-allocating) version of `metropolis_minmax`, filling existing arrays

Run a Metropolis sampler to estimate the extrema of a finite-range source
distribution `dist` using samples drawn from that distribution -- e.g.,
estimate zircon saturation and eruption ages from a distribution of zircon
crystallization ages.

### Examples
```julia
metropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, 2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
```
"""
function metropolis_minmax!(tmindist::DenseArray, tmaxdist::DenseArray, lldist::DenseArray, acceptancedist::BitVector, dist::Collection{<:Number}, mu::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}; burnin::Integer=0)
    @assert eachindex(tmindist) == eachindex(tmaxdist) == eachindex(lldist) == eachindex(acceptancedist) 

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9

    # Initial step sigma for Gaussian proposal distributions
    youngest, oldest = minimum(mu), maximum(mu)
    Δt = (oldest - youngest) + sqrt(sigma[argmin(mu)]^2 + sigma[argmax(mu)]^2)
    tminstep = tmaxstep = Δt / length(mu)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = youngest - sigma[argmin(mu)]
    tmaxₚ = tmax = oldest + sigma[argmax(mu)]

    # Log likelihood of initial proposal
    llₚ = ll = dist_ll(dist, mu, sigma, tmin, tmax) + prior_ll(mu, sigma, tmin, tmax)

    # Burnin
    for i=1:burnin
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tminstep*randn())
        (r > 0.5) && (tminₚ += tmaxstep*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu, sigma, tminₚ, tmaxₚ) + prior_ll(mu, sigma, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ-tmax)*stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist, tmaxdist, lldist, acceptancedist)
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tminstep*randn())
        (r > 0.5) && (tminₚ += tmaxstep*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu, sigma, tminₚ, tmaxₚ) + prior_ll(mu, sigma, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ-tmax)*stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            acceptancedist[i]=true
        end
        tmindist[i] = tmin
        tmaxdist[i] = tmax
        lldist[i] = ll
    end
    return tmindist, tmaxdist, lldist, acceptancedist
end
function metropolis_minmax!(tmindist::DenseArray{T}, tmaxdist::DenseArray{T}, tlldist::DenseArray{T}, lldist::DenseArray{T}, acceptancedist::BitVector, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; burnin::Integer=0, tllprior=Uniform(0,minimum(value.(age68.(analyses)))), lossprior=Uniform(0,100), method=:projection) where {T <: AbstractFloat}
    @assert eachindex(tmindist) == eachindex(tmaxdist) == eachindex(tlldist) == eachindex(lldist) == eachindex(acceptancedist)
    @assert (method === :projection || method === :bivariate) "Allowed methods are `:projection` or `:bivariate`"

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9

    # Process input analyses
    ellipses = Ellipse.(analyses)
    ages68 = log.(one(T) .+ (ellipses .|> e->e.y₀))./value(λ238U)
    ages = @. upperintercept(zero(T), ellipses)

    # Initial step sigma for Gaussian proposal distributions
    youngest, oldest = minimum(ages), maximum(ages)
    Δt = (oldest.val - youngest.val) + sqrt(+ oldest.err^2 + youngest.err^2)
    tminstep = tmaxstep = Δt / length(analyses)
    tllstep = youngest.val/50

    # Use oldest and youngest upper intercept ages for initial proposal
    tminₚ = tmin = value(youngest)
    tmaxₚ = tmax = value(oldest)
    tllₚ = tll = zero(T)

    # Ensure priors for time and amount of lead loss are appropriately truncated
    tllprior = truncated(tllprior, 0, minimum(ages68))
    lossprior = truncated(lossprior, 0, 100) # percent

    # Log likelihood of initial proposal
    ll = logpdf(tllprior, tll)
    if method === :projection
        ll += dist_ll(dist, ages, tmin, tmax) 
    elseif method === :bivariate
        ll += dist_ll(dist, analyses, tmin, tmax, tll) 
    end
    ll += prior_ll(ages, tmin, tmax)
    for i in eachindex(ages, ages68)
        loss = 100*max(one(T) - (ages68[i] - tll) / (value(ages[i]) - tll), zero(T))
        ll += logpdf(lossprior, loss)
    end
    llₚ = ll

    # Burnin
    for i = 1:burnin
        tminₚ, tmaxₚ, tllₚ = tmin, tmax, tll
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tminstep * randn()
        elseif r < 0.70
            tmaxₚ += tmaxstep * randn()
        else
            tllₚ += tllstep * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(tllₚ, ellipses)
        llₚ = logpdf(tllprior, tllₚ)
        if method === :projection
            llₚ += dist_ll(dist, ages, tminₚ, tmaxₚ)
        elseif method === :bivariate
            llₚ += dist_ll(dist, analyses, tminₚ, tmaxₚ, tllₚ)
        end
        llₚ += prior_ll(ages, tminₚ, tmaxₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - tllₚ) / (value(ages[i]) - tllₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            tll = tllₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist, tlldist)
        tminₚ, tmaxₚ, tllₚ = tmin, tmax, tll
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tminstep * randn()
        elseif r < 0.70
            tmaxₚ += tmaxstep * randn()
        else
            tllₚ += tllstep * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(tllₚ, ellipses)
        llₚ = logpdf(tllprior, tllₚ)
        if method === :projection
            llₚ += dist_ll(dist, ages, tminₚ, tmaxₚ)
        elseif method === :bivariate
            llₚ += dist_ll(dist, analyses, tminₚ, tmaxₚ, tllₚ)
        end
        llₚ += prior_ll(ages, tminₚ, tmaxₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - tllₚ) / (value(ages[i]) - tllₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tminstep = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmaxstep = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            tll = tllₚ
            acceptancedist[i]=true
        end
        tmindist[i] = tmin
        tmaxdist[i] = tmax
        tlldist[i] = tll
        lldist[i] = ll
    end
    return tmindist, tmaxdist, tlldist, lldist, acceptancedist
end
