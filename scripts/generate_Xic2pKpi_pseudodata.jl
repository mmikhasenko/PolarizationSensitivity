using Plots
using Plots.StatsBase
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK
using DelimitedFiles

using PolarizationSensitivity

import PolarizationSensitivity: intensity, interference, ellh
intensity(σs; pars) = intensity(σs, isobars; pars=pars)

#Set-up model
tbs = tbs_Ξc2pKπ
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv,Λ1520_pc, Λ1520_pv)

#Dalitz plot for genpars
let pars = [1.3, 1 - 1im, 1.2 - 0.5im, 2 + 0.3im, 2 - 0.6im, 2 + 1im]
    σ3v = range(lims3(tbs.ms)..., length=100)
    σ1v = range(lims1(tbs.ms)..., length=100)
    cal = [
        (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
        NaN :
        intensity(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
    heatmap(σ1v, σ3v, cal, colorbar=false,
        color=cgrad(:viridis, scale=:exp),
        xlab="σ₁ ≡ m²(Kπ) (GeV)",
        ylab="σ₃ ≡ m²(pK) (GeV)")
        savefig("Dalitz.pdf")
end


const genpars = [1.3, 1 - 1im, 1.2 - 0.5im, 2 + 0.3im, 2 - 0.6im, 2 + 1im]


function generate_sample(Intensity; Nev)
    s = flatDalitzPlotSample(tbs.ms; Nev=Nev)
    ws = Intensity.(s)
    maxweight = max(ws...)
    check(w) = w/maxweight > rand()
    s[check.(ws)]
end

@time generate_sample(σs->intensity(σs; pars=genpars); Nev=20_000)

const gs_large = generate_sample(σs->intensity(σs; pars=genpars); Nev=100_000);

writedlm(joinpath("data","sims","sample_Kstar=1.3,1-1im_Lambda=1.2-0.5im, 2+0.3im_Delta=2-0.6im,2+1im.txt"),
    [getproperty.(gs_large,:σ1) getproperty.(gs_large,:σ3)])

#plot Dalitz of data
histogram2d(getproperty.(gs_large,:σ1), getproperty.(gs_large,:σ3), bins=50)

