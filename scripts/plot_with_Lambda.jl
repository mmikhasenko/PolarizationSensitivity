using Plots
using PartialWaveFunctions
using Parameters

theme(:wong)


let (two_j,two_τ,two_λ) = (1,1,-1), cosθ = 0.3
    PartialWaveFunctions.wignerd_doublearg(two_j,two_τ,two_λ,cosθ), -sin(acos(cosθ)/2)
end

A(σ; m,Γ) = 1/(m^2-σ-1im*m*Γ)

δ(ν,λ) = (ν==λ ? 1 : 0)

H2(two_λ) = (two_λ==1 ? 1.0 : -1.0)
H1(two_τ) = (two_τ==1 ? 1.0 :  1.0)

# Λ1520 assume that s = 1/2
O(σ1,cosθ,two_ν,two_λ) = sum(
        δ(two_ν,two_λ)*H1(two_τ) * A(σ1; m=1.52, Γ=0.017) *
        wignerd_doublearg(1,two_τ,two_λ,cosθ)*H2(two_λ)
            for two_τ in [-1,1])

Intensity(σ1,cosθ) = sum(abs2(O(σ1,cosθ,two_ν,two_λ))
    for two_ν in [-1,1], two_λ in [-1,1])

let
    σ1v = 1.0:0.01:2.8
    cosθv = -1:0.01:1
    calv = [Intensity(σ1,cosθ) for cosθ in cosθv, σ1 in σ1v]
    heatmap(σ1v, cosθv, calv, title="intensity",
        xlab="m²(pK) (GeV)", ylab="cosθ")
end

struct test
    x
    y
    z
end

# 
@with_kw struct ThreeBodySystem
    m0::Float64
    m1::Float64
    m2::Float64
    m3::Float64
end
sumsq(tbs::ThreeBodySystem) = tbs.m0^2+tbs.m1^2+tbs.m2^2+tbs.m3^2
# 
# 
@with_kw struct Invariants
    σ1::Float64
    σ2::Float64
    σ3::Float64
end
function Invariants(tbs::ThreeBodySystem;σ1=-1.0,σ2=-1.0,σ3=-1.0)
    sign(σ1)+sign(σ2)+sign(σ3)!=1 && error("the method works with TWO invariants given")
    σ3 < 0 && return Invariants(σ1=σ1,σ2=σ2,σ3=sumsq(tbs)-σ1-σ2)
    σ1 < 0 && return Invariants(σ2=σ2,σ3=σ3,σ1=sumsq(tbs)-σ2-σ3)
    σ2 < 0 && return Invariants(σ3=σ3,σ1=σ1,σ2=sumsq(tbs)-σ3-σ1)
end
#
# 
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function cosθ23(σs::Invariants,tbs::ThreeBodySystem)
    i=2; j=3; k=1;
    m2s = [tbs.m1^2, tbs.m2^2, tbs.m3^2]
    s = tbs.m0^2
    σsv = [σs.σ1, σs.σ2, σs.σ3]
    # 
    EE4σ = (σsv[k]+m2s[i]-m2s[j])*(s-σsv[k]-m2s[k])
    pp4σ = sqrt(λ(σsv[k],m2s[i],m2s[j])*λ(s,σsv[k],m2s[k]))
    rest = σsv[j]-m2s[k]-m2s[i]
    return (2σsv[k]*rest-EE4σ)/pp4σ
end

tbs = ThreeBodySystem(m0=2.46867, m1=0.13957, m2=0.49367, m3=0.938)

Intensity_σσ(σ1,σ3) =
    Intensity(σ1,cosθ23(Invariants(tbs,σ1=σ1,σ3=σ3),tbs))

Invariants(tbs;σ1=1.5^2,σ3=1.0)

# Ξc(0)→π(1)K(2)p(3)
# σ1 = (p2+p3)^2 : σ1 > (m2+m3)^2
# σ1 = (p0-p1)^2 : σ1 < (m0-m1)^2
Kibble(σs::Invariants,tbs::ThreeBodySystem) = λ(λ(tbs.m0^2,tbs.m1^2,σs.σ1),
                                                λ(tbs.m0^2,tbs.m2^2,σs.σ2),
                                                λ(tbs.m0^2,tbs.m3^2,σs.σ3))

let
    σ1v = range((tbs.m2+tbs.m3)^2, (tbs.m0-tbs.m1)^2, length=110)
    σ3v = range((tbs.m1+tbs.m2)^2, (tbs.m0-tbs.m3)^2, length=100)
    #
    calv = [(Kibble(Invariants(tbs,σ1=σ1,σ3=σ3),tbs) > 0 ?
        NaN :
        Intensity_σσ(σ1,σ3)) for σ3 in σ3v, σ1 in σ1v]
    #
    plot(layout = grid(2,2, widths=(0.3,0.7), heights=(0.7,0.3)), size=(800,800))
    heatmap!(sp=2,σ1v, σ3v, calv, title="intensity",
        xlab="m²(pK) (GeV)", ylab="m²(Kπ) (GeV)", colorbar=false)
    calv_noNaN = [isnan(c) ? 0.0 : c for c in calv]
    plot!(sp=4, σ1v, vcat(sum(calv_noNaN, dims=1)...), link=:x)
    plot!(sp=1, vcat(sum(calv_noNaN, dims=2)...), σ3v, link=:y)
end

using QuadGK

# integration
plojection_σ1(σ1) = quadgk(σ3->
    (Kibble(Invariants(tbs,σ1=σ1,σ3=σ3),tbs) > 0 ?
        0 :
        Intensity_σσ(σ1,σ3)),
    (tbs.m1+tbs.m2)^2, (tbs.m0-tbs.m3)^2)[1]

plojection_σ3(σ3) = quadgk(σ1->
    (Kibble(Invariants(tbs,σ1=σ1,σ3=σ3),tbs) > 0 ?
        0 :
        Intensity_σσ(σ1,σ3)),
    (tbs.m2+tbs.m3)^2, (tbs.m0-tbs.m1)^2)[1]

plot(σ1->plojection_σ1(σ1), (tbs.m2+tbs.m3)^2, (tbs.m0-tbs.m1)^2)
plot(σ3->plojection_σ3(σ3), (tbs.m1+tbs.m2)^2, (tbs.m0-tbs.m3)^2)