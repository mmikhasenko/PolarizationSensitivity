
const tbs_Ξc2pKπ = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));

#
ρ1(σ1; ms) = sqrt(
    λ(σ1, ms.m0^2, ms.m1^2)*
    λ(σ1, ms.m2^2, ms.m3^2))/σ1
