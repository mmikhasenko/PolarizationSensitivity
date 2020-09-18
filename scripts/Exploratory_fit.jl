using Plots
using Optim

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles

# 
tbs = tbs_Ξc2pKπ
# 
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv,Λ1520_pc, Λ1520_pv)
Np = length(isobars)

import PolarizationSensitivity: intensity, interference, ellh, fit_data!
intensity(σs; pars) = intensity(σs, isobars; pars=pars)
interference(σs; i,j) = interference(σs, isobars; i=i, j=j)
ellh(pars;data,H) = ellh(pars,isobars;data=data,H=H)


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

# reading
data = let Nreduced = 300
M = readdlm(joinpath("data","sims","sample_Kstar=1.3,1-1im_Lambda=1.2-0.5im, 2+0.3im_Delta=2-0.6im,2+1im.txt"))
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end

settings = Dict(
    "H_matrix" => H,
    "data"=> data,
    "Natt"=>1,
    "show_trace"=>false)

#Calculate normalized genpars
const genpars′ = genpars./sqrt(μ(genpars; H=H)/length(data))
@assert μ(genpars′; H=H) ≈ length(data)

#Perform fit
@time fit_data!(settings);

#Save fit settings to file
using JLD2
@save joinpath("data","fit_Np=6_Natt=10_pseudodata.jld2") settings


