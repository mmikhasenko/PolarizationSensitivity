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

# converged = filter(x->x.st, tfr);

# print(converged)

# histogram(converged.min, bins=100)

# bestpars = converged[findmin(converged.min)[2]].pars
# print(bestpars)

#const LMs = NamedTuple{(:L,:M)}.([(1,1),(2,1),(4,1)])
#Np = length(LMs)

