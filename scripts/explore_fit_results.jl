using Plots
using Optim
using Statistics

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
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv,Λ1520_pc, Λ1520_pv)
const parity_map = [1,0,1,0,1,0]


Np = length(isobars)

import PolarizationSensitivity: intensity, interference, ellh
intensity(σs; pars) = intensity(σs, isobars; pars=pars)
interference(σs; i,j) = interference(σs, isobars; i=i, j=j)
ellh(pars;data,H) = ellh(pars,isobars;data=data,H=H)

#Load fit settings
@load joinpath("data","fits_Np=6_Natt=10_pseudodata.jld2") settings

H = settings["H_matrix"]

Optim.minimizer(settings["fit_results"][1])

tfr = Table(
    [(st = Optim.converged(fr), min = Optim.minimum(fr), pars = Optim.minimizer(fr))
        for fr in settings["fit_results"]])
print(tfr)
print(tfr.pars[1])


#1.) Sanity check
histogram(μ.(tfr.pars; H = H), bins = range(90,110,length = 20), title = "Histogram of μ (sanity check)")
savefig("muhisto.pdf")

#2.) Fraction of every isobar 
let Np = length(isobars)
    titles = ["K*_pc","K*_pv","Δ**_pc","Δ**_pv", "Λ**_pc", "Λ**_pv"]
    plot( title = "Fractions of isobars")
    for i in 1:6

        bins = range(0,100,length = 50)

        calv = real(getindex.(contributions.(tfr.pars; H=H),i,i))

        histogram!(calv, lab =  "$(titles[i]) mean = $(round(mean(calv); digits = 2))",
            bins = bins, seriescolor = i)
        exp_cal = real(getindex(contributions(genpars′; H=H),i,i))
        vline!([exp_cal], lab="$(round(exp_cal; digits =2)) expected", seriescolor = i)
    end
    savefig("Fractions_isobars.pdf")
end


#Draw Dalitz with each isobar for 1 fit attempt
let pars = [0 + 0im , 0.0 , 0.0, 0.0,0.0,0.0]
    l = length(pars)
    plot(layout=grid(3,3,heights=(0.5,0.5,0.5)), size=(1200,1200))
    for j in 1:l
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
        pars =  [0 + 0im , 0.0 , 0.0, 0.0,0.0,0.0]
    end
    savefig("Dalitz_plot.pdf")
end

#3.) Relative phase between K*_pc and K*_pv
histogram(relative_phase.(tfr.pars,1,2), 
bins = range(-2,2, length = 50), title = "Relative phase between K*_pc and K*_pv")
vline!([ relative_phase(genpars′,1,2) ], lab="expected")
savefig("Relative_phase.pdf")


#4.) Intensity in the PV and PC sectors
μ(genpars′; H=H)
let 
    bins = range(0,100, length = 50)
    plot(title = "PV/PC fractions")
    histogram!(real(μ_pc.(tfr.pars; H=H)),
        bins = bins, lab = "Fraction of PC", seriescolor = 1)
    vline!([ real(μ_pc(genpars′; H=H)) ], lab="expected PC", seriescolor = 1)

    histogram!(real(μ_pv.(tfr.pars; H=H)),
        bins = bins, lab = "Fraction of PV", seriescolor = 2)
    vline!([ real(μ_pv(genpars′; H=H)) ], lab="expected PV", seriescolor = 2)
    #savefig("Fraction_of_PCPV.pdf")
end

#5.) Averaged αs
const sample = flatDalitzPlotSample(tbs.ms; Nev = 10000)

@time histogram(α_avg.(I1_assym,tfr.pars;isobars=isobars, sample=sample),
 bins = range(-1,1, length = 50), title = "Average α(K*)")
 vline!([ α_avg(I1_assym,genpars′; isobars=isobars, sample=sample) ], lab="expected")
 savefig("AlphaK*_avg.pdf")

@time histogram(α_avg.(I2_assym, tfr.pars; isobars = isobars, sample=sample), 
bins = range(-1,1, length = 50), title = "Average α(Δ**)")
vline!([ α_avg(I2_assym,genpars′; isobars=isobars, sample=sample) ], lab="expected")
savefig("AlphaDelta_avg.pdf")

@time histogram(α_avg.(I2_assym, tfr.pars; isobars = isobars, sample=sample),
bins = range(-1,1, length = 50), title = "Average α(Λ*)")
vline!([ α_avg(I3_assym, genpars′; isobars=isobars, sample=sample) ], lab="expected")
savefig("AlphaLambda_avg.pdf")


α_avg(I1_assym, genpars .* [1,1,0,0,0,0] ; isobars = isobars, sample = sample)
α_avg(I2_assym, genpars .* [0,0,1,1,0,0] ; isobars = isobars, sample = sample)
α_avg(I3_assym, genpars .* [0,0,0,0,1,1] ; isobars = isobars, sample = sample)



