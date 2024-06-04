
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
function dist_ll(dist::Collection, mu::Collection, sigma::Collection, tmin::Number, tmax::Number)
    tmax >= tmin || return NaN
    any(isnan, mu) && return NaN
    any(x->!(x>0), sigma) && return NaN
    @assert issorted(mu)
    mu₋, sigma₋ = first(mu), first(sigma)
    mu₊, sigma₊ = last(mu), last(sigma)
    nbins = length(dist) - 1
    dt = abs(tmax-tmin)

    # Cycle through each datum in dataset
    loglikelihood = zero(float(eltype(dist)))
    @inbounds for j in eachindex(mu, sigma)
        μⱼ, σⱼ = mu[j], sigma[j]

        # Find equivalent index position of μⱼ in the `dist` array
        ix = (μⱼ - tmin) / dt * nbins + 1
        # If possible, prevent aliasing problems by interpolation
        if (σⱼ < dt / nbins) && 1 < ix < length(dist)
            # Interpolate corresponding distribution value
            f = floor(Int, ix)
            δ = ix - f
            likelihood = (dist[f+1]*δ + dist[f]*(1-δ)) / dt
        else
            # Otherwise, sum contributions from Gaussians at each point in distribution
            𝑖 = 1:length(dist)
            likelihood = zero(float(eltype(dist)))
            normconst = 1/(length(dist) * σⱼ * sqrt(2 * pi))
            @turbo for i in eachindex(dist, 𝑖)
                distx = tmin + dt * (𝑖[i] - 1) / nbins # time-position of distribution point
                # Likelihood curve follows a Gaussian PDF. Note: dt cancels
                likelihood += dist[i] * normconst * exp(-(distx - μⱼ)^2 / (2 * σⱼ * σⱼ))
            end
        end
        loglikelihood += log(likelihood)
    end
    # Calculate a weighted mean and examine our MSWD
    (wm, wsigma, mswd) = wmean(mu, sigma, corrected=true) # Height of MSWD distribution relative to height at MSWD = 1
    # (see Wendt and Carl, 1991, Chemical geology)
    f = length(mu) - 1
    Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1)) * (f > 0)
    Zf = max(min(Zf, 1.0), 0.0)
    @assert 0 <= Zf <= 1
    # To prevent instability / runaway of the MCMC for small datasets (low N),
    # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
    # the youngest-zircon interpretation at low Zf (MSWD far from one). The
    # penalty factors used here are determined by training against synthetic datasets.
    # In other words, these are just context-dependent prior distributions on tmax and tmin
    loglikelihood -= (2/log(1+length(mu))) * (          # Scaling factor that decreases with log number of data points (i.e., no penalty at high N)
      log((abs(tmin - wm)+wsigma)/wsigma)*Zf +          # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmax - wm)+wsigma)/wsigma)*Zf +          # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmin - mu₋)+sigma₋)/sigma₋)*(1-Zf) +     # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
      log((abs(tmax - mu₊)+sigma₊)/sigma₊)*(1-Zf) )     # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    return loglikelihood
end
function dist_ll(dist::Collection, analyses::Collection{<:Measurement}, tmin::Number, tmax::Number)
    tmax >= tmin || return NaN
    any(isnan, analyses) && return NaN
    any(x->!(err(x) > 0), analyses) && return NaN
    old = maximum(analyses)
    yng = minimum(analyses)
    nbins = length(dist) - 1
    dt = abs(tmax - tmin)

    # Cycle through each datum in dataset
    loglikelihood = zero(float(eltype(dist)))
    @inbounds for j in eachindex(analyses)
        dⱼ = analyses[j]
        μⱼ, σⱼ = val(dⱼ), err(dⱼ)

        # Find equivalent index position of μⱼ in the `dist` array
        ix = (μⱼ - tmin) / dt * nbins + 1
        # If possible, prevent aliasing problems by interpolation
        if (σⱼ < dt / nbins) && 1 < ix < length(dist)
            # Interpolate corresponding distribution value
            f = floor(Int, ix)
            δ = ix - f
            likelihood = (dist[f+1]*δ + dist[f]*(1-δ)) / dt
        else
            # Otherwise, sum contributions from Gaussians at each point in distribution
            𝑖 = 1:length(dist)
            likelihood = zero(float(eltype(dist)))
            normconst = 1/(length(dist) * σⱼ * sqrt(2 * pi))
            @turbo for i in eachindex(dist, 𝑖)
                distx = tmin + dt * (𝑖[i] - 1) / nbins # time-position of distribution point
                # Likelihood curve follows a Gaussian PDF. Note: dt cancels
                likelihood += dist[i] * normconst * exp(-(distx - μⱼ)^2 / (2 * σⱼ * σⱼ))
            end
        end
        loglikelihood += log(likelihood)
    end
    # Calculate a weighted mean and examine our MSWD
    (wm, mswd) = wmean(analyses, corrected=true)
    @assert wm.err > 0
    # Height of MSWD distribution relative to height at MSWD = 1
    # (see Wendt and Carl, 1991, Chemical geology)
    f = length(analyses) - 1
    Zf = exp((f / 2 - 1) * log(mswd) - f / 2 * (mswd - 1)) * (f > 0)
    Zf = max(min(Zf, 1.0), 0.0)
    @assert 0 <= Zf <= 1

    # To prevent instability / runaway of the MCMC for small datasets (low N),
    # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
    # the youngest-zircon interpretation at low Zf (MSWD far from one). The
    # penalties used here were determined by training against synthetic datasets.
    # In other words, these are just context-dependent prior distributions on tmax and tmin
    loglikelihood -= (2 / log(1 + length(analyses))) * (                    # Scaling factor that decreases with log number of analyses points (i.e., no penalty at high N)
      log((abs(tmin - wm.val) + wm.err) / wm.err) * Zf +            # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmax - wm.val) + wm.err) / wm.err) * Zf +            # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmin - yng.val) + yng.err) / yng.err) * (1 - Zf) +   # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
      log((abs(tmax - old.val) + old.err) / old.err) * (1 - Zf))    # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    return loglikelihood
end


"""
```julia
metropolis_min(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; burnin::Integer=0, t0prior=Uniform(0,minimum(age68.(analyses))), lossprior=Uniform(0,100))
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
metropolis_min(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; kwargs...) = metropolis_min(nsteps, dist, val.(data), err.(data); kwargs...)
function metropolis_min(nsteps::Integer, dist::Collection, mu::Collection, sigma::Collection; kwargs...)
    # Allocate ouput array
    tmindist = Array{float(eltype(mu))}(undef,nsteps)
    # Run Metropolis sampler
    return metropolis_min!(tmindist, nsteps, dist, mu, sigma; kwargs...)
end
function metropolis_min(nsteps::Integer, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; kwargs...) where {T}
    # Allocate ouput arrays
    tmindist = Array{T}(undef,nsteps)
    t0dist = Array{T}(undef,nsteps)
    # Run Metropolis sampler
    metropolis_min!(tmindist, t0dist, nsteps, dist, analyses; kwargs...)
    return tmindist, t0dist
end


"""
```julia
metropolis_min!(tmindist::DenseArray, nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)
metropolis_min!(tmindist::DenseArray, t0dist::DenseArray, nsteps::Integer, dist::Collection, analyses::Collection{<:UPbAnalysis}; burnin::Integer=0) where {T}
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
function metropolis_min!(tmindist::DenseArray, nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)
    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9
    # Sort the dataset from youngest to oldest
    sI = sortperm(mu)
    mu_sorted = mu[sI] # Sort means
    sigma_sorted = sigma[sI] # Sort uncertainty
    youngest, oldest = first(mu_sorted), last(mu_sorted)

    # Step sigma for Gaussian proposal distributions
    dt = oldest - youngest + first(sigma_sorted) + last(sigma_sorted)
    tmin_step = dt / length(mu)
    tmax_step = dt / length(mu)
    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = youngest - first(sigma_sorted)
    tmaxₚ = tmax = oldest + last(sigma_sorted)

    # Log likelihood of initial proposal
    llₚ = ll = dist_ll(dist, mu_sorted, sigma_sorted, tmin, tmax)

    # Burnin
    for i=1:burnin
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tmin_step*randn())
        (r > 0.5) && (tminₚ += tmax_step*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu_sorted, sigma_sorted, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ-tmax)*stepfactor
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
        (r < 0.5) && (tmaxₚ += tmin_step*randn())
        (r > 0.5) && (tminₚ += tmax_step*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu_sorted, sigma_sorted, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ-tmax)*stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
        end
        tmindist[i] = tmin
    end
    return tmindist
end
function metropolis_min!(tmindist::DenseArray{T}, t0dist::DenseArray{T}, nsteps::Integer, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; burnin::Integer=0, t0prior=Uniform(0,minimum(age68.(analyses))), lossprior=Uniform(0,100)) where {T}
    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9
    # Sort the dataset from youngest to oldest

    # These quantities will be used more than once
    t0ₚ = t0 = 0.0
    ellipses = Ellipse.(analyses)
    ages68 = log.(one(T) .+ (ellipses .|> e->e.y₀))./val(λ238U)
    ages = similar(ellipses, Measurement{T})
    @. ages = upperintercept(t0ₚ, ellipses)
    youngest = minimum(ages)
    oldest = maximum(ages)
    t0step = youngest.val/50
    # t0prior = Uniform(0, youngest.val)
    t0prior = truncated(t0prior, 0, minimum(age68.(analyses)))
    lossprior = truncated(lossprior, 0, 100)

    # Initial step sigma for Gaussian proposal distributions
    dt = sqrt((oldest.val - youngest.val)^2 + oldest.err^2 + youngest.err^2)
    tmin_step = tmax_step = dt / length(analyses)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = val(youngest)
    tmaxₚ = tmax = val(oldest)

    # Log likelihood of initial proposal
    ll = dist_ll(dist, ages, tmin, tmax) + logpdf(t0prior, t0)
    for i in eachindex(ages, ages68)
        loss = 100*max(one(T) - (ages68[i] - t0) / (val(ages[i]) - t0), zero(T))
        ll += logpdf(lossprior, loss)
    end
    llₚ = ll

    # Burnin
    for i = 1:burnin
        tminₚ, tmaxₚ, t0ₚ = tmin, tmax, t0
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tmin_step * randn()
        elseif r < 0.70
            tmaxₚ += tmax_step * randn()
        else
            t0ₚ += t0step * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(t0ₚ, ellipses)
        llₚ = dist_ll(dist, ages, tminₚ, tmaxₚ)
        llₚ += logpdf(t0prior, t0ₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - t0ₚ) / (val(ages[i]) - t0ₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            t0 = t0ₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist, t0dist)
        tminₚ, tmaxₚ, t0ₚ = tmin, tmax, t0
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tmin_step * randn()
        elseif r < 0.70
            tmaxₚ += tmax_step * randn()
        else
            t0ₚ += t0step * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(t0ₚ, ellipses)
        llₚ = dist_ll(dist, ages, tminₚ, tmaxₚ)
        llₚ += logpdf(t0prior, t0ₚ)
        for i in eachindex(ages, ages68)
            loss = 100*max(one(T) - (ages68[i] - t0ₚ) / (val(ages[i]) - t0ₚ), zero(T))
            llₚ += logpdf(lossprior, loss)
        end
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            t0 = t0ₚ
        end
        tmindist[i] = tmin
        t0dist[i] = t0
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
metropolis_minmax(nsteps::Integer, dist::Collection, data::Collection{<:Measurement}; kwargs...) = metropolis_minmax(nsteps, dist, val.(data), err.(data); kwargs...)
function metropolis_minmax(nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; kwargs...)
    # Allocate ouput arrays
    acceptancedist = falses(nsteps)
    lldist = Array{float(eltype(dist))}(undef,nsteps)
    tmaxdist = Array{float(eltype(mu))}(undef,nsteps)
    tmindist = Array{float(eltype(mu))}(undef,nsteps)
    # Run metropolis sampler
    return metropolis_minmax!(tmindist, tmaxdist, lldist, acceptancedist, nsteps, dist, mu, sigma; kwargs...)
end
function metropolis_minmax(nsteps::Integer, dist::Collection{T}, analyses::Collection{<:UPbAnalysis{T}}; kwargs...) where T
    # Allocate ouput arrays
    acceptancedist = falses(nsteps)
    lldist = Array{T}(undef,nsteps)
    t0dist = Array{T}(undef,nsteps)
    tmaxdist = Array{T}(undef,nsteps)
    tmindist = Array{T}(undef,nsteps)
    # Run metropolis sampler
    return metropolis_minmax!(tmindist, tmaxdist, t0dist, lldist, acceptancedist, nsteps, dist, analyses; kwargs...)
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
function metropolis_minmax!(tmindist::DenseArray, tmaxdist::DenseArray, lldist::DenseArray, acceptancedist::BitVector, nsteps::Integer, dist::Collection, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0)
    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9
    # Sort the dataset from youngest to oldest
    sI = sortperm(mu)
    mu_sorted = mu[sI] # Sort means
    sigma_sorted = sigma[sI] # Sort uncertainty
    youngest, oldest = first(mu_sorted), last(mu_sorted)

    # Step sigma for Gaussian proposal distributions
    dt = oldest - youngest + first(sigma_sorted) + last(sigma_sorted)
    tmin_step = tmax_step = dt / length(mu)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = youngest - first(sigma_sorted)
    tmaxₚ = tmax = oldest + last(sigma_sorted)

    # Log likelihood of initial proposal
    llₚ = ll = dist_ll(dist, mu_sorted, sigma_sorted, tmin, tmax)

    # Burnin
    for i=1:nsteps
        # Adjust upper or lower bounds
        tminₚ, tmaxₚ = tmin, tmax
        r = rand()
        (r < 0.5) && (tmaxₚ += tmin_step*randn())
        (r > 0.5) && (tminₚ += tmax_step*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu_sorted, sigma_sorted, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ-tmax)*stepfactor
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
        (r < 0.5) && (tmaxₚ += tmin_step*randn())
        (r > 0.5) && (tminₚ += tmax_step*randn())
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        llₚ = dist_ll(dist, mu_sorted, sigma_sorted, tminₚ, tmaxₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ-tmin)*stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ-tmax)*stepfactor
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
function metropolis_minmax!(tmindist::DenseArray{T}, tmaxdist::DenseArray{T}, t0dist::DenseArray{T}, lldist::DenseArray{T}, acceptancedist::BitVector, nsteps::Integer, dist::Collection{T}, analyses::Collection{UPbAnalysis{T}}; burnin::Integer=0) where {T}
    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9
    # Sort the dataset from youngest to oldest

    # These quantities will be used more than once
    t0ₚ = t0 = 0.0
    ellipses = Ellipse.(analyses)
    ages = similar(ellipses, Measurement{T})
    @. ages = upperintercept(t0ₚ, ellipses)
    youngest = minimum(ages)
    oldest = maximum(ages)
    t0step = youngest.val/50
    t0prior = Uniform(0, youngest.val)

    # Initial step sigma for Gaussian proposal distributions
    dt = sqrt((oldest.val - youngest.val)^2 + oldest.err^2 + youngest.err^2)
    tmin_step = tmax_step = dt / length(analyses)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = val(youngest)
    tmaxₚ = tmax = val(oldest)

    # Log likelihood of initial proposal
    ll = llₚ = dist_ll(dist, ages, tmin, tmax) + logpdf(t0prior, t0)

    # Burnin
    for i = 1:burnin
        tminₚ, tmaxₚ, t0ₚ = tmin, tmax, t0
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tmin_step * randn()
        elseif r < 0.70
            tmaxₚ += tmax_step * randn()
        else
            t0ₚ += t0step * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(t0ₚ, ellipses)
        llₚ = dist_ll(dist, ages, tminₚ, tmaxₚ)
        llₚ += logpdf(t0prior, t0ₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            t0 = t0ₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i in eachindex(tmindist, t0dist)
        tminₚ, tmaxₚ, t0ₚ = tmin, tmax, t0
        # Adjust upper or lower bounds, or Pb-loss time
        r = rand()
        if r < 0.35
            tminₚ += tmin_step * randn()
        elseif r < 0.70
            tmaxₚ += tmax_step * randn()
        else
            t0ₚ += t0step * randn()
        end
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        @. ages = upperintercept(t0ₚ, ellipses)
        llₚ = dist_ll(dist, ages, tminₚ, tmaxₚ)
        llₚ += logpdf(t0prior, t0ₚ)
        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ - ll)
            if tminₚ != tmin
                tmin_step = abs(tminₚ - tmin) * stepfactor
            end
            if tmaxₚ != tmax
                tmax_step = abs(tmaxₚ - tmax) * stepfactor
            end

            ll = llₚ
            tmin = tminₚ
            tmax = tmaxₚ
            t0 = t0ₚ
            acceptancedist[i]=true
        end
        tmindist[i] = tmin
        tmaxdist[i] = tmax
        t0dist[i] = t0
        lldist[i] = ll
    end
    return tmindist, tmaxdist, t0dist, lldist, acceptancedist
end
