using Isoplot, Plots, Distributions, DelimitedFiles

# Example Detrital zircon U-Pb dataset (Karlstrom et al. 2018)
cd(@__DIR__)
data = readdlm("data/sixtymileA.csv", ',')

# Turn into UPbAnalysis objects
analyses = UPbAnalysis.(eachcol(data)...,)

cc = plot(xlabel="²⁰⁷Pb/²³⁵U", ylabel="²⁰⁶Pb/²³⁸U", framestyle=:box)
plot!(cc, analyses, color=:black, alpha=0.05, label="")

# Screen for discordance
t = first.(age.(analyses)) .< 1000
@info "Excluding $(count(.!t)) Paleoproterozoic analyses"
analyses = analyses[t]

# Plot in Wetherill concordia space
plot!(cc, analyses, color=:darkblue, alpha=0.5, label="")
ylims!(cc, 0, last(ylims(cc)))
xlims!(cc, 0, last(xlims(cc)))
concordiacurve!(cc) # Add concordia curve
savefig(cc, "concordia.pdf")
display(cc)

## --- Bayesian Pb-loss-aware deposition age estimation
nsteps = 10^6
tmindist, t0dist = metropolis_min(nsteps, UniformDistribution, analyses, burnin=10^4,lossprior=Normal(0,30))
## ---
tpbloss = CI(t0dist)
tdepo = CI(tmindist)
display(tdepo)
display(tpbloss)

# Add to concordia plot
I = rand(1:length(tmindist), 200) # Pick 100 random samples from the posterior distribution
concordialine!(cc, t0dist[I], tmindist[I], truncate=true, color=:darkred, alpha=0.05, label="Deposition: $tdepo Ma", legend=:bottomright) # Add to Concordia plot
savefig("concordiadepo.pdf")
display(cc)
tdepo

## --- Histogram of distribution of deposition age

te = histogram(tmindist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Deposition: $tdepo Ma", color=:darkred, fill=true, alpha=0.75, linealpha=0.1, framestyle=:box)
ylims!(te, 0, last(ylims()))
savefig(te, "depositionage.pdf")
display(te)


## --- Histogram of distribution of time of Pb-loss

tpb = histogram(t0dist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Pb-loss: $tpbloss Ma", color=:darkblue, fill=true, alpha=0.65, linealpha=0.1, framestyle=:box)
xlims!(tpb, 0, last(xlims()))
ylims!(tpb, 0, last(ylims()))
savefig(tpb, "pbloss.pdf")
display(tpb)

## --- Plot stacked

h = plot(cc,te,tpb, layout=(3,1), size=(500,1000), left_margin=(8,:mm))
savefig(h, "depositionestimation.pdf")
display(h)

## --- Standard deposition age estimation

nsteps = 10^6
tmindist68 = metropolis_min(nsteps, UniformDistribution, first.(age.(analyses)), burnin=10^4)
tdep68 = CI(tmindist68)
display(tdep68)

## ---
