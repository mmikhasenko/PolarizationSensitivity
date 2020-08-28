using Plots
using PartialWaveFunctions
theme(:wong)

# Λ(3/2+) → p K
A(cosθ,two_ν,two_λ; L) = 
    clebschgordan_doublearg(2*L,0,3,two_ν,1,two_ν)*
        wignerd_doublearg(3,two_ν,two_λ,cosθ)*
        clebschgordan_doublearg(2,0,1,two_λ,3,two_λ)

I(cosθ; L) = sum(abs2, A(cosθ,two_ν,two_λ; L=L) for two_ν in [-1,1], two_λ in [-1,1])

[clebschgordan_doublearg(2,0,3,two_ν,1,two_ν) for two_ν in [-1,1]]
[clebschgordan_doublearg(4,0,3,two_ν,1,two_ν) for two_ν in [-1,1]]

let
    plot()
    cosθv = range(-1,1,length=100)
    calv = I.(cosθv; L=1)
    plot!(cosθv, calv/sum(calv))
    calv = I.(cosθv; L=2)
    plot!(cosθv, calv/sum(calv))
end

let 
    plot()
    for (i,(two_ν,two_λ)) in enumerate(Iterators.product([-1,1],[-1,1]))
        plot!(cosθ->A(cosθ,two_ν,two_λ; L=1)*sqrt(3), -1, 1, seriescolor=i, lab="($two_ν,$two_λ)")
        plot!(cosθ->A(cosθ,two_ν,two_λ; L=2)*sqrt(5), -1, 1, seriescolor=i, ls=:dash, lw=2, lab="")
    end
    plot!()
end
atan(0,-1)