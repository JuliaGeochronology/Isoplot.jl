# Renne et al. 2011 reply 10.1016/j.gca.2011.06.021
const κFCs = 1.6417E-03 ± 0.0045E-03
const λ40Kϵ = 0.5757E-4 ± 0.0017E-4 # 1/Myr
const λ40Kβ = 4.9548E-4 ± 0.0134E-4 # 1/Myr
const λ40K = λ40Kβ + λ40Kϵ
export λ40K, λ40Kβ, λ40Kϵ, κFCs

# σκσλϵ =  7.1903E-13
# σκσλβ = -6.5839E-13
# σλϵσλβ = -3.4711E-14 # 1/Myr^2
#
# λ40KΣ = [ stdev(κFCs)^2  σκσλϵ  σκσλβ
#           σκσλϵ stdev(λ40Kϵ)^2 σλϵσλβ
#           σκσλβ  σλϵσλβ  stdev(λ40Kβ)^2]
