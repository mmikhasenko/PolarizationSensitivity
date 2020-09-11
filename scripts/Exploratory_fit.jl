using Plots
using Optim

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles

# 
tbs = tbs_Ξc2pKπ
# 
isobars = (Kst872_pc, Kst872_pv, Λ1520_pc, Δ1232_pc)
Np = length(isobars)

import PolarizationSensitivity: intensity, interference, ellh, fit_data!
intensity(σs; pars) = intensity(σs, isobars; pars=pars)
interference(σs; i,j) = interference(σs, isobars; i=i, j=j)
ellh(pars;data,H) = ellh(pars,isobars;data=data,H=H)

#Calculate integrals
Φ0 = quadgk(σ1->ρ1(σ1; tbs.ms), lims1(tbs.ms)...)[1]

#Calculate all matrix elements Interf_ij
const s0 = flatDalitzPlotSample(tbs.ms; Nev = 10000)
H = Matrix{Complex{Float64}}(undef,Np,Np)
@time for i in 1:Np
    for j in 1:Np
        H[i,j] = (Φ0/length(s0))*sum(interference.(s0; i=i,j=j))
    end
end

plot(heatmap(log.(abs.(real.(H)))),
    heatmap(imag.(H)))

# reading
data = let Nreduced = 100
    M = readdlm(joinpath("data","sims","sample_Kstar=1,1_Lambda=1.1_Delta=3.9.txt"))
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end

# LIKELIHOOD FUNCTION 
const genpars = [1, 1, 1.1, 3.9]
const genpars′ = genpars./sqrt(μ(genpars; H=H)/length(data))
@assert μ(genpars′; H=H) ≈ length(data)


let Np = 4
    plot()
    for j in 1:4
        ranges = sqrt.(length(data) ./ [H[i,i] for i in 1:4])
        ph = Vector{Float64}()
        e = Vector{Float64}()
        gen_pars = [1 + 0im, 1, 1.1, 3.9]
        for i in 1:100
            phase = (i-1)/99*2π - π
            modified_pars = copy(gen_pars)
            modified_pars[j] *= cis(phase)
            append!(ph,phase)
            append!(e,ellh(modified_pars; data=data, H=H))
        end
        plot!(ph,e, lab = "$j")
    end
    plot!()
end


settings = Dict(
    "H_matrix" => H,
    "data"=> data,
    "Natt"=>10,
    "show_trace"=>false)


@time fit_data!(settings);

using JLD2
@save joinpath("data","fits3_Np=4_Natt=10_pseudodata.jld2") settings

# Optim.minimizer(settings["fit_results"][1])

# tfr = Table(
#     [(st = Optim.converged(fr), min = Optim.minimum(fr), pars = Optim.minimizer(fr))
#         for fr in settings["fit_results"]])
# print(tfr)

# converged = filter(x->x.st, tfr);

# print(converged)

# histogram(converged.min, bins=100)

# bestpars = converged[findmin(converged.min)[2]].pars
# print(bestpars)

#const LMs = NamedTuple{(:L,:M)}.([(1,1),(2,1),(4,1)])
#Np = length(LMs)

