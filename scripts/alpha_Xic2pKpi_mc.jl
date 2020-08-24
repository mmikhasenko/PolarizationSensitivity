using Plots
using Plots.StatsBase
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK
theme(:wong)

const tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
#
#                   p    K    π    Ξc
const tbs_parities_pc = ['+', '-', '-', '+'];
const tbs_parities_pv = ['+', '-', '-', '-'];

# chain-1: K*(872)
lsh(s,σ) = BW(σ, 0.89176, 0.05)
const Kst872_pc = decay_chain(1, lsh;
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

rp = randomPoint(tbs)
O(rp.σs, rp.two_λs; pars = [0.1,0.1,0.3,0.4]) # O_λ^ν(σ1,σ2)
# 
I(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars) for two_λ in [-1,1], two_ν in [-1,1])
#

let pars = [1, 1, 1, 1]
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

# 
# 
# plotting α1
# 
# 

I1(σs, two_ν; pars) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars) for two_λ in [-1,1])

I2(σs, two_ν; pars) = sum(abs2,
    sum(
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat12(σs, tbs.ms^2)) * O(σs,
            ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
            pars=pars) for two_ν′ in [-1,1])
        for two_λ in [-1,1])
# 
I3(σs, two_ν; pars) = sum(abs2,
    sum((two_ν==two_ν′ ? -1 : 1)*
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat31(σs, tbs.ms^2)) * O(σs,
            ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
            pars=pars) for two_ν′ in [-1,1])
        for two_λ in [-1,1])

let pars=[1, 1, 1.1, 3.9im]
    I2(rp.σs, 1; pars=pars) + I2(rp.σs, -1; pars=pars) ≈
        I1(rp.σs, 1; pars=pars) + I1(rp.σs, -1; pars=pars) ≈
        I3(rp.σs, 1; pars=pars) + I3(rp.σs, -1; pars=pars)
end

const s = flatDalitzPlotSample(tbs.ms, Nev=50_000);
mid(x) = (x[1:end-1]+x[2:end]) / 2

const bestpars = [1, 1.0, 1, 1];

p1 = let pars = bestpars
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600), link=:x, lab="")
    # 
    σ1v = getproperty.(s,:σ1)
    calv⁺ = I1.(s,  1; pars=pars)
    calv⁻ = I1.(s, -1; pars=pars)
    stephist!(sp=1, σ1v, weights = calv⁺ + calv⁻, ylabel = "Intensity", bins=70, lab="")
    # 
    h⁺ = fit(Histogram, σ1v, weights(calv⁺), nbins=70)
    h⁻ = fit(Histogram, σ1v, weights(calv⁻), nbins=70)
    αv = (h⁺.weights - h⁻.weights) ./ (h⁺.weights + h⁻.weights)
    plot!(sp=2, mid(h⁺.edges[1]), αv, ylims=(-1,1), ylab="α₁",
            xlabel = "m²(Kπ)", lab="")
    # length(h⁺.edges[1])
end

p2 = let pars = bestpars
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600), link=:x)
    # 
    σ3v = getproperty.(s,:σ3)
    calv⁺ = I3.(s,  1; pars=pars)
    calv⁻ = I3.(s, -1; pars=pars)
    stephist!(sp=1, σ3v, weights = calv⁺ + calv⁻, ylabel = "Intensity", bins=70, lab="")
    # 
    h⁺ = fit(Histogram, σ3v, weights(calv⁺), nbins=70)
    h⁻ = fit(Histogram, σ3v, weights(calv⁻), nbins=70)
    αv = (h⁺.weights - h⁻.weights) ./ (h⁺.weights + h⁻.weights)
    plot!(sp=2, mid(h⁺.edges[1]), αv, ylims=(-1,1), ylab="α₂",
            xlabel = "m²(pK)", lab="")
    # length(h⁺.edges[1])
end

p3 = let pars = bestpars
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600), link=:x)
    # 
    σ2v = getproperty.(s,:σ2)
    calv⁺ = I2.(s,  1; pars=pars)
    calv⁻ = I2.(s, -1; pars=pars)
    stephist!(sp=1, σ2v, weights = calv⁺ + calv⁻, ylabel = "Intensity", bins=70, lab="")
    # 
    h⁺ = fit(Histogram, σ2v, weights(calv⁺), nbins=70)
    h⁻ = fit(Histogram, σ2v, weights(calv⁻), nbins=70)
    αv = (h⁺.weights - h⁻.weights) ./ (h⁺.weights + h⁻.weights)
    plot!(sp=2, mid(h⁺.edges[1]), αv, ylims=(-1,1), ylab="α₃",
            xlabel = "m²(πp)", lab="")
    # length(h⁺.edges[1])
end

plot(p1,p2,p3, size=(1000,600), layout=grid(1,3))
savefig(joinpath("plots", "three_alphas_1.0_1.0_1.0_1.0.pdf"))