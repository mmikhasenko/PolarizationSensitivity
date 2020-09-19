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
const datafile = joinpath("data","sims","sample_Kstar=1.3,1-1im_Delta=1.2-0.5im,2+0.3im_Lambda=2-0.6im,2+1im.txt")
const data = let Nreduced = 100
    M = readdlm(datafile)
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end
const genpars = parse_values_from_datafile_name(datafile, Np)
const genpars′ = genpars./sqrt(μ(genpars; H=H)/length(data))

const H = JLD2.read(JLD2.jldopen(joinpath("data", "fit_Np=6_Natt=10_pseudodata.jld2")), "settings")["H_matrix"]

function profile_llh(pars)
    pars′ = pars./sqrt(μ(pars; H=H)/length(data))
    ellh(pars′, isobars; data=data, H=H)
end

let
    plot()
    ϕv = range(-π,π,length=50)
    for p in 1:Np
        mask = [p==i for i in 1:Np]
        calv = [profile_llh(genpars′.*(iszero.(mask) + mask.*cis(ϕ))) for ϕ in ϕv]
        plot!(ϕv, calv, xlab="ϕ of PV sector", ylab="likelihood", lab="par #$(p)")
    end
    plot!()
end
savefig(joinpath("plots", "profiled_llh_vs_pars_phase.pdf"))
