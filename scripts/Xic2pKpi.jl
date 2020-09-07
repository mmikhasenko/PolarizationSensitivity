using Plots
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK


tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
#
print(tbs)
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
Kst872_pv.two_ls
Kst872_pc.two_ls
Kst872_pv.two_LS
Kst872_pc.two_LS
print(Kst872_pc)
print(Kst872_pv)


# chain-2: Delta**
Δ1232 = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs,
    parity = '+', Ps = tbs_parities_pc)

# chain-3: Lambda**
Λ1520  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
# Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)

O1(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv,Λ1520,Δ1232)))

#print(tbs.ms^2)

function O2(σs,two_λs; pars) 
    two_j = 1
    two_ν = two_λs[4]
    two_λ = two_λs[1]
    
    x = sum(wignerd_doublearg(two_j,two_ν,two_ν′,cosθhat12(σs, tbs.ms^2))*
            sum(c*amplitude(σs,ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),dc)
                for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv,Λ1520,Δ1232)))
            for two_ν′ in [-1,1])
    return x
end

function O3(σs,two_λs; pars) 
    two_j = 1
    two_ν = two_λs[4]
    two_λ = two_λs[1]
    
    x = sum((-1)^(two_ν - two_ν′)*wignerd_doublearg(two_j,two_ν,two_ν′,cosθhat31(σs, tbs.ms^2))*
            sum(c*amplitude(σs,ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),dc)
                for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv,Λ1520,Δ1232)))
            for two_ν′ in [-1,1])
    return x
end

rp = randomPoint(tbs) 
O1(rp.σs, rp.two_λs; pars = [0.1,0.1,0.3,0.4]) # O_λ^ν(σ1,σ2)
O2(rp.σs, rp.two_λs; pars = [0.1,0.1,0.3,0.4]) # O_λ^ν(σ1,σ2)
O3(rp.σs, rp.two_λs; pars = [0.1,0.1,0.3,0.4]) # O_λ^ν(σ1,σ2)
#print(rp.σs)
#print(typeof(rp.two_λs))
# 
I1(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2(O1(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars)) for two_λ in [-1,1], two_ν in [-1,1])

I2(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2(O2(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars)) for two_λ in [-1,1], two_ν in [-1,1])

I3(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2(O3(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars)) for two_λ in [-1,1], two_ν in [-1,1])
#
#print(tbs.ms)
print(ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0)[1])

let 
    σ3v = range(lims3(tbs.ms)..., length=100)
    #σ3v = range(lims3(tbs.ms)..., length=100)
    σ1v = range(lims1(tbs.ms)..., length=100)
    cal = [
        (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms) > 0 ?
        NaN :
        I1(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=[1, 0, 1, 2])) for σ3 in σ3v, σ1 in σ1v]
    heatmap(σ1v, σ3v, cal, colorbar=false,
        color=cgrad(:viridis, scale=:exp),
        xlab="σ₁ ≡ m²(Kπ) (GeV)",
        ylab="σ₃ ≡ m²(pK) (GeV)")
#savefig("Dalitz_plot.png")
end

# 
# plotting α1
# 
#
integrate_over_σ3_as_function_of_σ1(σ1; integrand) = quadgk(σ3->
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
        0 :
        integrand(Invariants(tbs.ms,σ1=σ1,σ3=σ3))),
    lims3(tbs.ms)...)[1]

integrate_over_σ3_as_function_of_σ2(σ2; integrand) = quadgk(σ3->
    (Kibble(Invariants(tbs.ms,σ2=σ2,σ3=σ3),tbs.ms^2) > 0 ?
        0 :
        integrand(Invariants(tbs.ms,σ2=σ2,σ3=σ3))),
    lims3(tbs.ms)...)[1]

integrate_over_σ1_as_function_of_σ3(σ3; integrand) = quadgk(σ1->
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
        0 :
        integrand(Invariants(tbs.ms,σ3=σ3,σ1=σ1))),
    lims1(tbs.ms)...)[1]

#quadgk(x->sin(x),(-1,1)...)

I1_assym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O1(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
   pars=pars)) for two_λ in [-1,1])

I2_assym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O2(σs,
   ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
  pars=pars)) for two_λ in [-1,1])

I3_assym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O3(σs,
   ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
  pars=pars)) for two_λ in [-1,1])

I1_assym(rp.σs, 1;pars=[1, 1, 1.1, 3.9im])
I2_assym(rp.σs, 1;pars=[1, 1, 1.1, 3.9im])
I3_assym(rp.σs, 1;pars=[1, 1, 1.1, 3.9im]) 

I1(rp.σs;pars=[1, 1, 1.1, 3.9im])
I2(rp.σs;pars=[1, 1, 1.1, 3.9im])
I3(rp.σs;pars=[1, 1, 1.1, 3.9im])


function α1(σ1; pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1_assym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1_assym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α2(σ2;pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ3_as_function_of_σ2(σ2; integrand = σs->I2_assym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ3_as_function_of_σ2(σ2; integrand = σs->I2_assym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α3(σ3;pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ1_as_function_of_σ3(σ3; integrand = σs->I3_assym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ1_as_function_of_σ3(σ3; integrand = σs->I3_assym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end


α1(rp.σs.σ1; pars=[1, 1, 0, 0])
α2(rp.σs.σ2; pars=[1, 1, 1.1, 3.9im])
α3(rp.σs.σ3; pars=[1, 1, 1.1, 3.9im])

let pars = [1, 1, 1, 1]
    σ1v = range(lims1(tbs.ms)..., length = 100)
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600))
    # 
    calv = integrate_over_σ3_as_function_of_σ1.(σ1v; integrand = σs->I1(σs; pars=pars))
    plot!(sp=1, σ1v, calv, legend = false, title = "K* → Kπ resonance",
            xlabel = "m²(Kπ)", ylabel = "Intensity")
    # 
    calv = α1.(σ1v; pars=pars)
    plot!(sp=2, σ1v, calv, ylims=(-1,1),
            xlabel = "m²(Kπ)", ylabel = "α(K*)")
    calv = α1.(σ1v; pars=[1,1,0,0])
    plot!(sp = 2, σ1v, calv,  ylims=(-1,1), legend = false)
#    savefig("alfa_Kstar.png")
end

let pars = [1, 1, 1, 1]
    σ2v = range(lims2(tbs.ms)..., length = 100)
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600))
    # 
    calv = integrate_over_σ3_as_function_of_σ2.(σ2v; integrand = σs->I2(σs; pars=pars))
    plot!(sp=1, σ2v, calv,legend = false, title = "Δ** → pπ resonance",
            xlabel = "m²(pπ)", ylabel = "Intensity")
    # 
    calv = α2.(σ2v; pars=pars)
    plot!(sp=2, σ2v, calv, ylims=(-1,1),legend = false,
            xlabel = "m²(pπ)", ylabel = "α(Δ**)")
    calv = α2.(σ2v; pars=[0,0,1,0])
    plot!(sp=2, σ2v, calv, ylims=(-1,1), legend = false)
#    savefig("alfa_Deltastarstar.png")
end

let pars = [1, 1, 1, 1]
    σ3v = range(lims3(tbs.ms)..., length = 100)
    plot(layout=grid(2,1, heights=(0.7,0.3)), size=(500,600))
    # 
    calv = integrate_over_σ1_as_function_of_σ3.(σ3v; integrand = σs->I3(σs; pars=pars))
    plot!(sp=1, σ3v, calv,legend = false, title = "Λ* → pK resonance",
            xlabel = "m²(pK)", ylabel = "Intensity")
    # 
    calv = α3.(σ3v; pars=pars)
    plot!(sp=2, σ3v, calv, ylims=(-1,1),
            xlabel = "m²(pK)", ylabel = "α(Λ*)")
    calv = α3.(σ3v; pars=[0,0,0,1])
    plot!(sp=2, σ3v, calv, ylims=(-1,1), legend = false)
#    savefig("alfa_Lambdastar.png")
end



# θhat12 > 0
# θhat23 > 0
# θhat31 > 0
# 
# θhat21 clockwise < 0
# θhat32 clockwise < 0
# θhat13 clockwise < 0

cosθhat12(rp.σs, tbs.ms^2)
cosθhat23(rp.σs, tbs.ms^2)
cosθhat31(rp.σs, tbs.ms^2)


wignerd_doublearg(two_j,two_m1,two_m2,cosθhat12) # will work fine

wignerd_doublearg(two_j,two_m1,two_m2,cos(-acos(cosθhat21))) # incorrect

wignerd_doublearg(two_j,two_m1,two_m2,cosθhat21) = (-1)^{m1-m2} * 
            wignerd_doublearg(two_j,two_m1,two_m2,cosθhat12)

# α1 => chain 1 => d(cosθhat_1(3)) = (-1)^{...} d(cosθhat31)
# α2 => chain 2 => d(cosθhat_2(3)) = d(cosθhat23)