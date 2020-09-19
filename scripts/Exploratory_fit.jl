using Plots
using Optim

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles
using LinearAlgebra
using Statistics
#
#
#  
const tbs = tbs_Ξc2pKπ
#
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv,Λ1520_pc, Λ1520_pv)
const Np = length(isobars)
#
plot(tbs.ms, σs->interference(σs, isobars; i=1, j=1))
#
const datafile = joinpath("data","sims","sample_Kstar=1.3,1-1im_Delta=1.2-0.5im,2+0.3im_Lambda=2-0.6im,2+1im.txt")
const data = let Nreduced = 1000
    M = readdlm(datafile)
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end

const genpars = parse_values_from_datafile_name(datafile, Np)

#Calculate all matrix elements Interf_ij
const s0 = flatDalitzPlotSample(tbs.ms; Nev = 10000)
H = Matrix{Complex{Float64}}(undef,Np,Np)
@time for i in 1:Np
    for j in 1:Np
        H[i,j] = (Φ0/length(s0))*sum(interference.(s0; i=i,j=j))
    end
end

#Plot matrix
plot(heatmap(log.(abs.(real.(H)))),
    heatmap(imag.(H)))

settings = Dict(
    "H_matrix" => H,
    "data" => data,
    "Natt" => 1,
    "show_trace" => true,

#Calculate normalized genpars
const genpars′ = genpars./sqrt(μ(genpars; H=H)/length(data))
@assert μ(genpars′; H=H) ≈ length(data)

#Perform fit
@time fit_data!(settings);

#Save fit settings to file
using JLD2
@save joinpath("data","fit_Np=6_Natt=10_pseudodata.jld2") settings


