using Plots
using PartialWaveFunctions
using Parameters

theme(:wong)

A(σ; m,Γ) = 1/(m^2-σ-1im*m*Γ)

δ(ν,λ) = (ν==λ ? 1 : 0)

H2(two_λ) = (two_λ==1 ? 1.0 : -1.0)
H1(two_τ) = (two_τ==1 ? 1.0 :  1.0)

O(σ1,cosθ,two_ν,two_λ) = sum(
        δ(two_ν,two_τ)*H1(two_τ) * A(σ1; m=1.52, Γ=0.017) *
        wignerd_doublearg(1,two_τ,two_λ,cosθ)*H2(two_λ)
            for two_τ in [-1,1])
                
Intensity(σ1,cosθ) = sum(abs2(O(σ1,cosθ,two_ν,two_λ))
    for two_ν in [-1,1], two_λ in [-1,1])



α() = sum()