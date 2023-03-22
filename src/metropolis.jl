function dist_ll(
        dist::Collection,
        data::Collection{<:Measurement},
        tmin::Number,
        tmax::Number,
    )
    # Define some frequently used variables
    old = maximum(data)
    yng = minimum(data)
    nbins = length(dist) - 1
    dt = abs(tmax - tmin)

    # Cycle through each datum in dataset
    loglikelihood = zero(float(eltype(dist)))
    @inbounds for j in eachindex(data)
        dⱼ = data[j]
        μⱼ, σⱼ = val(dⱼ), err(dⱼ)

        # Find equivalent index position of μⱼ in the `dist` array
        ix = (μⱼ - tmin) / dt * nbins + 1
        # If possible, prevent aliasing problems by interpolation
        if (σⱼ < dt / nbins) && ix > 1 && ix < length(dist)
            # Interpolate corresponding distribution value
            f = floor(Int, ix)
            δ = ix - f
            likelihood = (dist[f+1]*δ + dist[f]*(1-δ)) / dt
            # Otherwise, sum contributions from Gaussians at each point in distribution
        else
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
    (wm, mswd) = awmean(data)
    # Height of MSWD distribution relative to height at MSWD = 1
    # (see Wendt and Carl, 1991, Chemical geology)
    f = length(data) - 1
    Zf = exp((f / 2 - 1) * log(mswd) - f / 2 * (mswd - 1)) * (f > 0)
    # To prevent instability / runaway of the MCMC for small datasets (low N),
    # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
    # the youngest-zircon interpretation at low Zf (MSWD far from one). The
    # penalties used here were determined by training against synthetic datasets.
    # In other words, these are just context-dependent prior distributions on tmax and tmin
    loglikelihood -= (2 / log(1 + length(data))) * (                    # Scaling factor that decreases with log number of data points (i.e., no penalty at high N)
      log((abs(tmin - wm.val) + wm.err) / wm.err) * Zf +            # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmax - wm.val) + wm.err) / wm.err) * Zf +            # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
      log((abs(tmin - yng.val) + yng.err) / yng.err) * (1 - Zf) +   # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
      log((abs(tmax - old.val) + old.err) / old.err) * (1 - Zf))    # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    return loglikelihood
end


function metropolis_min(nsteps::Int, dist::Collection, data::Collection{UPbAnalysis{T}}; burnin::Integer=0) where {T}
    # Allocate ouput arrays
    tmindist = Array{T}(undef,nsteps)
    t0dist = Array{T}(undef,nsteps)
    # Run Metropolis sampler
    metropolis_min!(tmindist, t0dist, nsteps, dist, data; burnin)
    return tmindist, t0dist
end

function metropolis_min!(
        tmindist::DenseArray,
        t0dist::DenseArray,
        nsteps::Integer,
        dist::Collection,
        data::Collection{<:UPbAnalysis};
        burnin::Integer = 0,
    )
    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
    stepfactor = 2.9
    # Sort the dataset from youngest to oldest

    # These quantities will be used more than once
    t0 = 0.0
    ellipses = ellipse.(data)
    ages = upperintercept.(t0, ellipses)
    youngest = minimum(ages)
    oldest = maximum(ages)
    t0step = youngest.val/50
    t0prior = Uniform(0, youngest.val)

    # Initial step sigma for Gaussian proposal distributions
    dt = sqrt((oldest.val - youngest.val)^2 + oldest.err^2 + youngest.err^2)
    tmin_step = tmax_step = dt / length(data)

    # Use oldest and youngest zircons for initial proposal
    tminₚ = tmin = youngest.val
    tmaxₚ = tmax = oldest.val

    # Log likelihood of initial proposal
    ll = llₚ = dist_ll(dist, ages, tmin, tmax) + logpdf(t0prior, t0)

    # Burnin
    for i = 1:burnin
        # Adjust upper and lower bounds
        tminₚ = tmin + tmin_step * randn()
        tmaxₚ = tmax + tmax_step * randn()
        t0ₚ = t0 + t0step * randn()
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
    @inbounds for i = 1:nsteps
        # Adjust upper and lower bounds
        tminₚ = tmin + tmin_step * randn()
        tmaxₚ = tmax + tmax_step * randn()
        t0ₚ = t0 + t0step * randn()
        # Flip bounds if reversed
        (tminₚ > tmaxₚ) && ((tminₚ, tmaxₚ) = (tmaxₚ, tminₚ))

        # Calculate log likelihood for new proposal
        for j in eachindex(ellipses, ages)
            ages[j] = upperintercept(t0ₚ, ellipses[j])
        end
        # @. ages = upperintercept(t0ₚ, ellipses)
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
        tmindist[i] = tmin
        t0dist[i] = t0

    end
    return tmindist
end
