using Isoplot, Plots

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
analyses = analyses[discordance.(analyses) .< 0.2]

# Plot in Wetherill concordia space
hdl = plot(xlabel="²⁰⁷Pb/²³⁵U", ylabel="²⁰⁶Pb/²³⁸U", framestyle=:box)
plot!(hdl, analyses, color=:darkblue, alpha=0.3, label = "")
concordiacurve!(hdl) # Add concordia curve
savefig(hdl, "concordia.svg")
display(hdl)

## --- Bayesian Pb-loss-aware eruption age estimation

nsteps = 100000
tmindist, t0dist = metropolis_min(nsteps, HalfNormalDistribution, analyses; burnin=10000)

using VectorizedStatistics
m = vmean(tmindist)
l = m - vpercentile(tmindist, 2.5)
u = vpercentile(tmindist, 97.5) - m
println("Eruption age: $m +$u/-$l Ma (95% CI)")


## --- Histogram of distribution of eruption age

h = histogram(tmindist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Eruption age", framestyle=:box)
ylims!(h, 0, ylims()[2])
savefig(h, "EruptionAge.svg")
display(h)


## --- Show eruption age relative to distribution of upper intercepts

uis = upperintercept.(vmean(t0dist), analyses)
h = rankorder(Isoplot.val.(uis), Isoplot.err.(uis))
plot!(h,1:length(uis),fill(m-l,length(uis)),fillto=m+u,color=:blue,fillalpha=0.5,linealpha=0, label="Model ($(round(m,digits=3)) +$(round(u,digits=3))/-$(round(l,digits=3)) Ma)")
plot!(h,1:length(uis),fill(m,length(uis)),linecolor=:black,linestyle=:dot,label="",legend=:topleft,fg_color_legend=:white,framestyle=:box)

# concordialine!(hdl, t₀, t₁, color=:darkred)


## --- Histogram of distribution of time of Pb-loss

h = histogram(t0dist, xlabel="Age [Ma]", ylabel="Probability Density", normalize=true, label="Time of Pb-loss", framestyle=:box)
xlims!(h, 0, xlims()[2])
ylims!(h, 0, ylims()[2])
savefig(h, "PbLoss.svg")
display(h)
