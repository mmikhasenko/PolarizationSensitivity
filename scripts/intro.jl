using DrWatson
@quickactivate "polarization_sensitivity"

using Plots
using ThreeBodyDecay

const mp = 0.938;
const mK = 0.49367;
const mπ = 0.13957;
const mΞc = 2.46867;

# create two-body system
tbs = ThreeBodySystem(mp,mK,mπ,mΞc;   # masses m1,m2,m3,m0
        two_jps=([ 1//2,   0,    0,  1//2] .|> x2,  # twice spin
                 [  '+',  '+',  '-',  '+']))        # parities

# chain-3, i.e. (1+2): Λs with the lowest ls, LS
Λ1520  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs);
Λ1690  = decay_chain(3, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs);
Λ1810  = decay_chain(3, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs);
#
# chain-1, i.e. (2+3): Pentaquarks with the lowest ls, LS
Ks892  = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05); two_s =   1|>x2, tbs=tbs);
#
Δ1232  = decay_chain(2, (s,σ)->BW(σ, 1.232,  0.112); two_s = 1/2|>x2, tbs=tbs);
#
A(σs,two_λs,cs) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(cs, (Λ1520, Λ1690, Λ1810, Ks892, Δ1232)))
I(σs,cs) = sum(abs2(A(σs,two_λs,cs)) for two_λs in itr(tbs.two_js))
#
I(σs)   = I(  σs,[1, 1.1, 3.9im, 2.2, 2.1im, -0.3im]) # set the couplings
#
dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
@show I(dpp.σs)

let
    #
    σ1v,σ2v,σ3v = flatDalitzPlotSample(tbs; Nev = 10000)
    #
    weights    = [I(  σs) for σs in zip(σ1v,σ2v,σ3v)]
    #
    histogram2d(σ1v,σ2v, weights = weights, bins=100)
end
