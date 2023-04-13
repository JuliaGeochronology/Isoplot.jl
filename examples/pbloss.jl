using Isoplot, Plots, VectorizedStatistics

# Example U-Pb dataset (MacLennan et al. 2020)
#       207/235  1σ abs   206/236     1σ abs     correlation
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

# Turn into UPbAnalysis objects
analyses = UPbAnalysis.(eachcol(data)...,)
# Screen for discordance
t = discordance.(analyses) .< 0.2
@info "Excluding $(count(.!t)) discordant analyses"
analyses = analyses[t]

# Plot in Wetherill concordia space
cc = plot(xlabel="²⁰⁷Pb/²³⁵U", ylabel="²⁰⁶Pb/²³⁸U", framestyle=:box)
plot!(cc, analyses, color=:darkblue, alpha=0.3, label="")
concordiacurve!(cc) # Add concordia curve
savefig(cc, "concordia.pdf")
display(cc)

## --- Log likelihood distribution / demonstration of burnin
nsteps = 10^6
tmindist, tmaxdist, t0dist, lldist, acceptancedist = metropolis_minmax(nsteps, UniformDistribution, analyses)
ll = plot(lldist, label="", xlabel="Step number", ylabel="Log likelihood", color=:darkblue, framestyle=:box, alpha=0.2)
savefig(ll, "lldist.pdf")
display(ll)

## --- Bayesian Pb-loss-aware eruption age estimation
nsteps = 10^7
tmindist, t0dist = metropolis_min(nsteps, UniformDistribution, analyses, burnin=10^4)
# ---
tpbloss = CI(t0dist)
terupt = CI(tmindist)
display(terupt)
display(tpbloss)

# Add to concordia plot
I = rand(1:length(tmindist), 200) # Pick 100 random samples from the posterior distribution
concordialine!(cc, t0dist[I], tmindist[I], color=:darkred, alpha=0.05, label="Eruption: $terupt Ma", legend=:bottomright) # Add to Concordia plot
savefig("concordiaterupt.pdf")
display(cc)
terupt

## --- Histogram of distribution of eruption age
te = histogram(tmindist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Eruption: $terupt Ma", color=:darkred, fill=true, alpha=0.75, linealpha=0.1, framestyle=:box)
ylims!(te, 0, last(ylims()))
savefig(te, "eruptionage.pdf")
display(te)


## --- Histogram of distribution of time of Pb-loss
tpb = histogram(t0dist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Pb-loss: $tpbloss Ma", color=:darkblue, fill=true, alpha=0.65, linealpha=0.1, framestyle=:box)
xlims!(tpb, 0, last(xlims()))
ylims!(tpb, 0, last(ylims()))
savefig(tpb, "pbloss.pdf")
display(tpb)

## --- Plot stacked
h = plot(cc,te,tpb, layout=(3,1), size=(500,1000), left_margin=(8,:mm))
savefig(h, "eruptionestimation.pdf")
display(h)


## --- Standard eruption age estimation

nsteps = 10^7
tmindist = metropolis_min(nsteps, UniformDistribution, first.(age.(analyses)), burnin=10^4)
teruptold = CI(tmindist)
display(teruptold)

## --- Stacked plot of screening options

# h = plot(ccfiltered,cc, layout=(2,1), size=(500,660), left_margin=(4,:mm))
# savefig(h, "filteredunfiltered.pdf")

## ---
