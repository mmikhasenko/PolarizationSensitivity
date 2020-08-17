using Plots
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK


tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
#
#                   p    K    π    Ξc
tbs_parities_pc = ['+', '-', '-', '+'];
tbs_parities_pv = ['+', '-', '-', '-'];

# chain-1: K*(872)
Kst872_pc = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Ξc -> K* p
# +  -> -  + # odd wave -> (-1)^L = - => L = 1
# K* J^P = 1^-  ⊗ p J^P = 1/2^+

Kst872_pv = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pv)
# Ξc -> K* p
# -  -> -  + # ever wave -> (-1)^L = - => L = 0
Kst872_pv.two_LS

# chain-2: Delta**
Δ1232 = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs,
    parity = '+', Ps = tbs_parities_pc)

# chain-3: Lambda**
Λ1520  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
# Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)

O(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv,Λ1520,Δ1232)))

rp = randomPoint(tbs)
O(rp.σs, rp.two_λs; pars = [0.1,0.1,0.3,0.4]) # O_λ^ν(σ1,σ2)
# 
I(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2(O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars)) for two_λ in [-1,1], two_ν in [-1,1])
#
let 
    σ3v = range(lims3(tbs.ms)..., length=100)
    σ1v = range(lims1(tbs.ms)..., length=100)
    cal = [
        (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms) > 0 ?
        NaN :
        I(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=[1, 1, 1.7, 0.3])) for σ3 in σ3v, σ1 in σ1v]
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

integrate_over_σ3_as_function_of_σ1(σ1; integrand) = quadgk(σ3->
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms) > 0 ?
        0 :
        integrand(Invariants(tbs.ms,σ1=σ1,σ3=σ3))),
    lims3(tbs.ms)...)[1]

I1(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars)) for two_λ in [-1,1])

function α1(σ1; pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1(σs, 1; pars=pars))
    O⁻ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

α1(rp.σs.σ1; pars=[1, 1, 1.1, 3.9im])

let pars = [1, 1, 1.7, 0.3]
    σ1v = range(lims1(tbs.ms)..., length=30)
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600))
    # 
    calv = integrate_over_σ3_as_function_of_σ1.(σ1v; integrand = σs->I(σs; pars=pars))
    plot!(sp=1, σ1v, calv,
            xlabel = "m²(Kπ)", ylabel = "Intensity")
    # 
    calv = α1.(σ1v; pars=pars)
    plot!(sp=2, σ1v, calv, ylims=(-1,1),
            xlabel = "m²(Kπ)", ylabel = "α")
end

# θhat12 > 0
# θhat23 > 0
# θhat31 > 0
# 
# θhat21 clockwise < 0
# θhat32 clockwise < 0
# θhat13 clockwise < 0

cosθhat12(rp.σs, tbs.ms^2)

wignerd_doublearg(two_j,two_m1,two_m2,cosθhat12) # will work fine

wignerd_doublearg(two_j,two_m1,two_m2,cos(-acos(cosθhat21))) # incorrect

wignerd_doublearg(two_j,two_m1,two_m2,cosθhat21) = (-1)^{m1-m2} * 
    wignerd_doublearg(two_j,two_m1,two_m2,cosθhat12)

# α1 => chain 1 => d(cosθhat_1(3)) = (-1)^{...} d(cosθhat31)
# α2 => chain 2 => d(cosθhat_2(3)) = d(cosθhat23)
