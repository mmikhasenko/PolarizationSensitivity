using Plots
using Optim

using TypedTables

using PolarizationSensitivity
using ThreeBodyDecay
using QuadGK
using DelimitedFiles
using PartialWaveFunctions

using JLD2

# 
tbs = tbs_Ξc2pKπ
# 
isobars = (Kst872_pc, Kst872_pv, Λ1520_pc, Δ1232_pc)

Np = length(isobars)

import PolarizationSensitivity: intensity, interference, ellh, I1_assym, α1_avg
intensity(σs; pars) = intensity(σs, isobars; pars=pars)
interference(σs; i,j) = interference(σs, isobars; i=i, j=j)
ellh(pars;data,H) = ellh(pars,isobars;data=data,H=H)
#I1_assym(σs, two_ν; pars) = I1_assym(σs, two_ν, isobars; pars = pars)
#α1_avg(pars; Φ0, sample) = α1_avg(pars, isobars; Φ0 = Φ0, sample = sample)

#Calculate integrals
Φ0 = quadgk(σ1->ρ1(σ1; tbs.ms), lims1(tbs.ms)...)[1]

#Calculate all matrix elements Interf_ij
const s0 = flatDalitzPlotSample(tbs.ms; Nev = 10000)
H = Matrix{Complex{Float64}}(undef,Np,Np)
@time for i in 1:Np
    for j in 1:Np
        H[i,j] = (Φ0/length(s0))*sum(interference.(s0; i=i,j=j))
    end
end


@load joinpath("data","fits4_Np=4_Natt=100_pseudodata.jld2") settings

Optim.minimizer(settings["fit_results"][1])

tfr = Table(
    [(st = Optim.converged(fr), min = Optim.minimum(fr), pars = Optim.minimizer(fr))
        for fr in settings["fit_results"]])
print(tfr)
print(tfr.pars[1])


#1.) Sanity check
histogram(μ.(tfr.pars; H = H), bins = range(90,110,length = 20), title = "Histogram of μ")
savefig("muhisto.pdf")

#2.) Fraction of every isobar 
histogram(real(getindex.(contributions.(tfr.pars; H=H),1,1)),
    bins = range(0,100,length = 50), title = "Fraction of K*_pc")
savefig("Fraction_of_K_pc.pdf")

histogram(real(getindex.(contributions.(tfr.pars; H=H),2,2)),
    bins = range(0,100,length = 50), title = "Fraction of K*_pv")
savefig("Fraction_of_K_pc.pdf")

histogram(real(getindex.(contributions.(tfr.pars; H=H),3,3)),
    bins = range(0,100,length = 50), title = "Fraction of Λ*")
savefig("Fraction_of_Lambda.pdf")

histogram(real(getindex.(contributions.(tfr.pars; H=H),4,4)),
    bins = range(0,100,length = 50), title = "Fraction of Δ**")
savefig("Fraction_of_Delta.pdf")


# Test μ calculated both ways
sum(contributions(tfr.pars[2]; H=H))
μ(tfr.pars[2]; H=H)

#Draw Dalitz with each isobar for 1 fit attempt
let pars = [0 + 0im , 0.0 , 0.0, 0.0]
    plot(layout=grid(2,2, heights=(0.5,0.5)), size=(1200,1200))
    for j in 1:4
        pars[j] = tfr.pars[1][j]
        label = μ(pars; H=H)
        σ3v = range(lims3(tbs.ms)..., length=100)
        σ1v = range(lims1(tbs.ms)..., length=100)
        cal = [
            (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
            NaN :
            intensity(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
        heatmap!(sp = j, σ1v, σ3v, cal, colorbar=true,
            color=cgrad(:viridis, scale=:exp),
            xlab="σ₁ ≡ m²(Kπ) (GeV)",
            ylab="σ₃ ≡ m²(pK) (GeV)", title = "μ = $label")
        pars =  [0 + 0im , 0.0 , 0.0, 0.0]
        savefig("Dalitz2.pdf")
    end
end

#3.) Relative phase between K* and Δ**
histogram(relative_phase.(tfr.pars,1,3), 
bins = range(-2,2, length = 50), title = "Relative phase between K* and Δ**")
savefig("Relative_phase.pdf")


#4.) Intensity in the PV and PC sectors
#Define μ_pc and μ_pv

histogram(real(μ_pc.(tfr.pars; H=H)),
bins = range(0,100, length = 50), title = "Fraction of PC")
savefig("Fraction_of_muPC.pdf")

histogram(real(μ_pv.(tfr.pars; H=H)),
bins = range(0,100, length = 50), title = "Fraction of PV")
savefig("Fraction_of_muPV.pdf")

#5.) Averaged αs
const sample = flatDalitzPlotSample(tbs.ms; Nev = 10000)
@time histogram(α1_avg.(tfr.pars;isobars=isobars, sample=sample),
 bins = range(-1,1, length = 50), title = "Average α(K*)")
 savefig("Alpha1_avg.pdf")

@time histogram(α2_avg.(tfr.pars; isobars = isobars, sample=sample), 
bins = range(-1,1, length = 50), title = "Average α(Λ*)")
savefig("Alpha2_avg.pdf")

@time histogram(α3_avg.(tfr.pars; isobars = isobars, sample=sample),
bins = range(-1,1, length = 50), title = "Average α(Δ**)")
savefig("Alpha3_avg.pdf")


I2_assym(σs, two_ν; pars, isobars) = sum(abs2,sum(
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat12(σs, tbs.ms^2)) * O(σs,
            ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
            pars=pars, isobars = isobars) for two_ν′ in [-1,1])
        for two_λ in [-1,1])

I3_assym(σs, two_ν; pars, isobars) = sum(abs2,
        sum((-1)^(two_ν - two_ν′)*
        wignerd_doublearg(1,two_ν,two_ν′,cosθhat31(σs, tbs.ms^2)) * O(σs,
        ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),
        pars=pars, isobars = isobars) for two_ν′ in [-1,1])
    for two_λ in [-1,1])


I1_assym(σs, two_ν; pars, isobars) = sum(abs2(O(σs,ThreeBodySpins(
        two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0);
        pars=pars, isobars = isobars)) for two_λ in [-1,1])

function α1_avg(pars; isobars, sample)
    N⁺ = (Φ0/length(sample))*sum(I1_assym.(sample, 1; pars=pars, isobars = isobars))
    N⁻ = (Φ0/length(sample))*sum(I1_assym.(sample, -1; pars=pars, isobars = isobars))
    return (N⁺ - N⁻) / (N⁺ + N⁻)
end
function α2_avg(pars; isobars,sample)
    N⁺ = (Φ0/length(sample))*sum(I2_assym.(sample, 1; pars=pars, isobars = isobars))
    N⁻ = (Φ0/length(sample))*sum(I2_assym.(sample, -1; pars=pars, isobars = isobars))
    return (N⁺ - N⁻) / (N⁺ + N⁻)
end

function α3_avg(pars; isobars,sample)
    N⁺ = (Φ0/length(sample))*sum(I3_assym.(sample, 1; pars=pars, isobars = isobars))
    N⁻ = (Φ0/length(sample))*sum(I3_assym.(sample, -1; pars=pars, isobars = isobars))
    return (N⁺ - N⁻) / (N⁺ + N⁻)
end


α1_avg(tfr.pars[1]; isobars=isobars, sample=sample)
α2_avg(tfr.pars[1]; isobars=isobars, sample=sample)
α3_avg(tfr.pars[1]; isobars=isobars, sample=sample)

# implement
import PolarizationSensitivity: α1
α1(pars) = α1(pars, isobars)
# α2(pars, isobars)
# α3(pars, isobars)

# add the plot for fit attemps
α1v = α1.(tfr.pars)
histogram(α1v, xlab="α₁")

# add tests
#a = -1-1im
#atan(imag(a) / real(a)) / π * 180 # ∈ [-π/2, π/2]
#atan(imag(a), real(a)) / π * 180 # ∈ [-π, π]
#arg(p[1]'*p[3]) # is phase of 3 - phase of 1

# 
# p3 = |p3|*exp(iϕ3)
# p1 = |p1|*exp(iϕ1)
# p1'*p3 = |p3|*|p1|*exp(iϕ3-iϕ1)

# converged = filter(x->x.st, tfr);

# print(converged)

# histogram(converged.min, bins=100)

# bestpars = converged[findmin(converged.min)[2]].pars
# print(bestpars)

#const LMs = NamedTuple{(:L,:M)}.([(1,1),(2,1),(4,1)])
#Np = length(LMs)

