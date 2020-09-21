using Optim
using JLD2

function load_results(hcycle;
    pathname = joinpath("data", "fit"),
    Natt = 50, Ndata = 1000, Np = 6)
return JLD2.read(
    JLD2.jldopen("$(pathname)_Ndata=$(Ndata)_Np=$(Np)_Natt=$(Natt)_hcycle=$(hcycle)_pseudodata.jld2"), "settings")["fit_results"]
end

allresults = vcat(load_results.(vcat(1:10,100:109))...);

let Natt = 50, Ndata = 1000, Np = 6, hcycle = 1
    @load joinpath("data","fit_Ndata=$(Ndata)_Np=$(Np)_Natt=$(Natt)_hcycle=$(hcycle)_pseudodata.jld2") settings
    settings["Natt"] = length(allresults)
    settings["fit_results"] = allresults
    # save
    Natt = length(allresults)
    @save joinpath("data","fit_Ndata=$(Ndata)_Np=$(Np)_Natt=$(Natt)_pseudodata.jld2") settings
end
