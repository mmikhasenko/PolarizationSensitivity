using Plots
using Plots.StatsBase
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK


const tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
#
#                   p    K    π    Ξc
const tbs_parities_pc = ['+', '-', '-', '+'];
const tbs_parities_pv = ['+', '-', '-', '-'];

# chain-1: K*(872)
const Kst872_pc = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Ξc -> K* p
# +  -> -  + # odd wave -> (-1)^L = - => L = 1
# K* J^P = 1^-  ⊗ p J^P = 1/2^+

const Kst872_pv = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pv)
# Ξc -> K* p
# -  -> -  + # ever wave -> (-1)^L = - => L = 0

# chain-2: Delta**
const Δ1232 = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs,
    parity = '+', Ps = tbs_parities_pc)

# chain-3: Lambda**
const Λ1520  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
# Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)

O(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv,Λ1520,Δ1232)))

# 
I(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars) for two_λ in [-1,1], two_ν in [-1,1])
#

let pars = [1, 1, 1.7, 2.3]
    σ3v = range(lims3(tbs.ms)..., length=100)
    σ1v = range(lims1(tbs.ms)..., length=100)
    cal = [
        (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
        NaN :
        I(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
    heatmap(σ1v, σ3v, cal, colorbar=false,
        color=cgrad(:viridis, scale=:exp),
        xlab="σ₁ ≡ m²(Kπ) (GeV)",
        ylab="σ₃ ≡ m²(pK) (GeV)")
end



const genpars = [1, 1, 1.1, 3.9im]

function generate_sample(Intensity; Nev)
    s = flatDalitzPlotSample(tbs.ms; Nev=Nev)
    ws = Intensity.(s)
    maxweight = max(ws...)
    check(w) = w/maxweight > rand()
    s[check.(ws)]
end

@time generate_sample(σs->I(σs; pars=genpars); Nev=20_000)


const gs_large = generate_sample(σs->I(σs; pars=genpars); Nev=100_000);

using DelimitedFiles
writedlm(joinpath("data","sims","sample_Kstar=1,1_Lambda=1.1_Delta=3.9im.txt"),
    [getproperty.(gs_large,:σ1) getproperty.(gs_large,:σ3)])

histogram2d(getproperty.(gs_large,:σ1), getproperty.(gs_large,:σ3), bins=50)

# reading
let
    M = readdlm(joinpath("data","sims","sample_Kstar=1,1_Lambda=1.1_Delta=3.9im.txt"))
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)]
end
