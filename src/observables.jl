# 
O(σs, two_λs; pars, isobars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, isobars))
#
interference(σs, isobars; i,j) = sum(
    conj(O(σs, two_λs; pars = delta(i; N=length(isobars)), isobars=isobars))*
         O(σs, two_λs; pars = delta(j; N=length(isobars)), isobars=isobars) 
    for two_λs in itr(tbs_Ξc2pKπ.two_js))
#
intensity(σs, isobars; pars) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars, isobars=isobars) for two_λ in [-1,1], two_ν in [-1,1])


#
const parity_map = [1,0,1,0,1,0]
μ_pc(pars;H) = μ(pars .* parity_map; H = H)
μ_pv(pars;H) = μ(pars .* iszero.(parity_map); H = H)


#Calculate integrals
Φ0 = quadgk(σ1->ρ1(σ1; tbs_Ξc2pKπ.ms), lims1(tbs_Ξc2pKπ.ms)...)[1]

I1_assym(σs, two_ν; pars, isobars) = sum(abs2(O(σs,ThreeBodySpins(
        two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0);
        pars=pars, isobars = isobars)) for two_λ in [-1,1])

I2_assym(σs, two_ν; pars, isobars) = sum(abs2,sum(
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat12(σs, tbs_Ξc2pKπ.ms^2)) * O(σs,
            ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
            pars=pars, isobars = isobars) for two_ν′ in [-1,1])
        for two_λ in [-1,1])

I3_assym(σs, two_ν; pars, isobars) = sum(abs2,
        sum((-1)^(two_ν - two_ν′)*
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat31(σs, tbs_Ξc2pKπ.ms^2)) * O(σs,
        ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
        pars=pars, isobars = isobars) for two_ν′ in [-1,1])
    for two_λ in [-1,1])

function α_avg(I_assym, pars; isobars, sample)
    N⁺ = (Φ0/length(sample))*sum(I_assym.(sample, 1; pars=pars, isobars = isobars))
    N⁻ = (Φ0/length(sample))*sum(I_assym.(sample, -1; pars=pars, isobars = isobars))
    return (N⁺ - N⁻) / (N⁺ + N⁻)
end
