using Plots
using Optim

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles
using LinearAlgebra
using Statistics
using JLD2
#
#
const tbs = tbs_Ξc2pKπ
#
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv, Λ1520_pc, Λ1520_pv)
const Np = length(isobars)
#
#
const datafile = joinpath("data","sims","sample_Kstar=1.3,1-1im_Delta=2-0.6im,2+1im_Lambda=1.2-0.5im,2+0.3im.txt")
const data = let Nreduced = 1000
    M = readdlm(datafile)
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end

# Calculate all matrix elements Interf_ij
# const H = integral_matrix_using_MC(isobars[1:1]; Nev=10_000)
const H = JLD2.read(JLD2.jldopen(joinpath("data", "integral_matrix_Np=6.jld2")), "settings")["H_matrix"]

const genpars = parse_values_from_datafile_name(datafile, Np)
const genpars′ = genpars./sqrt(μ(genpars; H=H)/length(data)) # normalized genpars
@assert μ(genpars′; H=H) ≈ length(data)

settings = Dict(
    "H_matrix" => H,
    "data" => data,
    "Natt" => 50,
    "show_trace" => false,
    "llh_tolerance" => 1e-4)

@time ellh(genpars′, isobars; data=data, H=H)

# # Perform fit
# @time fit_data!(settings, isobars);

# # Save fit settings to file
# let Natt = settings["Natt"], Ndata = length(data)
#     @save joinpath("data","fit_Ndata=$(Ndata)_Np=$(Np)_Natt=$(Natt)_pseudodata.jld2") settings
# end

# save every hour-cycle
for hcycle = 1:10
    println("starting $(hcycle) hcycle")
    # 
    fit_data!(settings, isobars);
    Natt = settings["Natt"]
    Ndata = length(data)
    #
    @save joinpath("data","fit_Ndata=$(Ndata)_Np=$(Np)_Natt=$(Natt)_hcycle=$(hcycle)_pseudodata.jld2") settings
end