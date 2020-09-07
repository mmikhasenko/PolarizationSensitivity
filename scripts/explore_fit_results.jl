using Plots
using Optim

using TypedTables

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles

using JLD2

@load joinpath("data","exp_pro","fits_Np=4_Natt=3_pseudodata.jld2") settings

Optim.minimizer(settings["fit_results"][1])

tfr = Table(
    [(st = Optim.converged(fr), min = Optim.minimum(fr), pars = Optim.minimizer(fr))
        for fr in settings["fit_results"]])
print(tfr)

# TODO:
# - move to utils
arg(z) = atan(imag(z), real(z))

# add tests
a = -1-1im
atan(imag(a) / real(a)) / π * 180 # ∈ [-π/2, π/2]
atan(imag(a), real(a)) / π * 180 # ∈ [-π, π]
arg(p[1]'*p[3]) # is phase of 3 - phase of 1

# implement
import PolarizationSensitivity: α1
α1(pars) = α1(pars, isobars)
# α2(pars, isobars)
# α3(pars, isobars)

# add the plot for fit attemps
α1v = α1.(tfr.pars)
histogram(α1v, xlab="α₁")

# 
# p3 = |p3|*exp(iϕ3)
# p1 = |p1|*exp(iϕ1)
# p1'*p3 = |p3|*|p1|*exp(iϕ3-iϕ1)


# converged = filter(x->x.st, tfr);

# print(converged)

# histogram(converged.min, bins=100)

# bestpars = converged[findmin(converged.min)[2]].pars
# print(bestpars)

#const LMs = NamedTuple{(:L,:M)}.([(1,1),(2,1),(4,1)])
#Np = length(LMs)

