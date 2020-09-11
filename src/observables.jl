# 
O(σs,two_λs; pars, isobars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, isobars))
#
interference(σs, isobars; i,j) = sum(
    conj(O(σs,two_λs; pars = delta(i;N=length(isobars)), isobars=isobars))*
         O(σs,two_λs;pars = delta(j;N=length(isobars)), isobars=isobars) 
    for two_λs in itr(tbs_Ξc2pKπ.two_js))
#
intensity(σs, isobars; pars) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars, isobars=isobars) for two_λ in [-1,1], two_ν in [-1,1])

I1_assym(σs, two_ν; pars) = sum(abs2(O(σs,ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0);
   pars=pars, isobars = isobars)) for two_λ in [-1,1])


function α1_avg(pars; sample)
    N⁺ = (Φ0/length(sample))*sum(I1_assym.(sample, 1; pars=pars))
    N⁻ = (Φ0/length(sample))*sum(I1_assym.(sample, -1; pars=pars))
    return (N⁺ - N⁻) / (N⁺ + N⁻)
end

