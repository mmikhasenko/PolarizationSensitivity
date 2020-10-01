using Plots
using PartialWaveFunctions
using ThreeBodyDecay
using QuadGK
using PolarizationSensitivity
using DelimitedFiles


const tbs = tbs_Ξc2pKπ
#                   p    K    π    Ξc
tbs_parities_pc = ['+', '-', '-', '+'];
tbs_parities_pv = ['+', '-', '-', '-'];

const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv, Λ1520_pc, Λ1520_pv)
const Np = length(isobars)
# II. load fitted data
# 
datafile = joinpath("data","sims","sample_Kstar=1.3,1-1im_Delta=2-0.6im,2+1im_Lambda=1.2-0.5im,2+0.3im.txt")
genpars = parse_values_from_datafile_name(datafile, Np)


O1(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, isobars))


function O2(σs,two_λs; pars) 
    two_j = 1
    two_ν = two_λs[4]
    two_λ = two_λs[1]
    
    x = sum(wignerd_doublearg(two_j,two_ν,two_ν′,cosθhat12(σs, tbs.ms^2))*
            sum(c*amplitude(σs,ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),dc)
                for (c, dc) in zip(pars, isobars))
            for two_ν′ in [-1,1])
    return x
end

function O3(σs,two_λs; pars) 
    two_j = 1
    two_ν = two_λs[4]
    two_λ = two_λs[1]
    
    x = sum((-1)^(two_ν - two_ν′)*wignerd_doublearg(two_j,two_ν,two_ν′,cosθhat31(σs, tbs.ms^2))*
            sum(c*amplitude(σs,ThreeBodySpins(two_h0=two_ν′, two_h1=two_λ, two_h2=0, two_h3=0),dc)
                for (c, dc) in zip(pars, isobars))
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
I1(σs; pars=[0.1, 0.1, 1.1, 3.9im, 5, 2]) = sum(abs2(O1(σs,
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
    σ1v = range(lims1(tbs.ms)..., length=100)
    cal = [
        (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
        NaN :
        I1(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=genpars)) for σ3 in σ3v, σ1 in σ1v]
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

I1_asym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O1(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
   pars=pars)) for two_λ in [-1,1])

I2_asym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O2(σs,
   ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
  pars=pars)) for two_λ in [-1,1])

I3_asym(σs, two_ν; pars=[1, 1, 1.1, 3.9im]) = sum(abs2(O3(σs,
   ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
  pars=pars)) for two_λ in [-1,1])

function α1(σ1; pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1_asym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ3_as_function_of_σ1(σ1; integrand = σs->I1_asym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α1_nonint(σs; pars)
    O⁺ = I1_asym(σs, 1; pars=pars)
    O⁻ = I1_asym(σs, -1; pars=pars)
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α2(σ2;pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ3_as_function_of_σ2(σ2; integrand = σs->I2_asym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ3_as_function_of_σ2(σ2; integrand = σs->I2_asym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α2_nonint(σs; pars)
    O⁺ = I2_asym(σs, 1; pars=pars)
    O⁻ = I2_asym(σs, -1; pars=pars)
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α3(σ3;pars=[1, 0.5, 1.1, 3.9im])
    O⁺ = integrate_over_σ1_as_function_of_σ3(σ3; integrand = σs->I3_asym(σs, 1; pars=pars))
    O⁻ = integrate_over_σ1_as_function_of_σ3(σ3; integrand = σs->I3_asym(σs, -1; pars=pars))
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

function α3_nonint(σs; pars)
    O⁺ = I3_asym(σs, 1; pars=pars)
    O⁻ = I3_asym(σs, -1; pars=pars)
    return (O⁺ - O⁻) / (O⁺ + O⁻)
end

let nbins = 100
    #Plot α1
    σ1v = range(lims1(tbs.ms)..., length = nbins)
    plot(layout=grid(2,3, heights=(0.7,0.3)), size=(1200,400), 
        bottom_margin=4Plots.PlotMeasures.mm)
    # 
    calv = integrate_over_σ3_as_function_of_σ1.(σ1v; integrand = σs->I1(σs; pars=genpars))
    plot!(sp=1, σ1v, calv, legend = false, title = "K* → Kπ resonance",
            xlabel = "m²(Kπ)", ylabel = "Intensity")
    # 
    calv = α1.(σ1v; pars=genpars)
    plot!(sp=4, σ1v, calv, ylims=(-1,1),
            xlabel = "m²(Kπ)", ylabel = "α(K*)")
    calv = α1.(σ1v; pars = ([1,1,0,0,0,0] .* genpars))
    plot!(sp = 4, σ1v, calv,  ylims=(-1,1), legend = false)
    #--------------------------------------------------------------------
    #Plot α2
    σ2v = range(lims2(tbs.ms)..., length = nbins)
    # 
    calv = integrate_over_σ3_as_function_of_σ2.(σ2v; integrand = σs->I2(σs; pars=genpars))
    plot!(sp= 2, σ2v, calv,legend = false, title = "Δ** → pπ resonance",
            xlabel = "m²(pπ)", ylabel = "Intensity")
    # 
    calv = α2.(σ2v; pars=genpars)
    plot!(sp=5, σ2v, calv, ylims=(-1,1),legend = false,
            xlabel = "m²(pπ)", ylabel = "α(Δ**)")
    calv = α2.(σ2v; (pars=[0,0,1,1,0,0] .* genpars) )
    plot!(sp=5, σ2v, calv, ylims=(-1,1), legend = false)
#    savefig("alfa_Kstar.png")
    #----------------------------------------------------------------------
    #Plot α3
    σ3v = range(lims3(tbs.ms)..., length = nbins)
    # 
    calv = integrate_over_σ1_as_function_of_σ3.(σ3v; integrand = σs->I3(σs; pars=genpars))
    plot!(sp=3, σ3v, calv,legend = false, title = "Λ* → pK resonance",
            xlabel = "m²(pK)", ylabel = "Intensity")
    # 
    calv = α3.(σ3v; pars=genpars)
    plot!(sp=6, σ3v, calv, ylims=(-1,1),
            xlabel = "m²(pK)", ylabel = "α(Λ*)")
    calv = α3.(σ3v; ( pars=[0,0,0,0,1,1] .* genpars) )
    plot!(sp=6, σ3v, calv, ylims=(-1,1), legend = false)
end
savefig("alphas.pdf")

let pars = [1.3, 1 - 1im, 2.0 - 0.6im, 2 + 1im, 1.2 - 0.5im, 2 + 0.3im]
    plot(layout = grid(1,3,widths=(0.29,0.29,0.42)), size=(3*400,350), bottom_margin=4Plots.PlotMeasures.mm)
    titles = ["α(K*)", "α(Δ**)", "α(Λ*)"]
    σ1v = range(lims1(tbs.ms)..., length = 100)
    σ3v = range(lims3(tbs.ms)..., length = 100)
   # cal = [
   # (α1_nonint(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
   cal = [
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
    NaN :
    α1_nonint(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]

    heatmap!(σ1v, σ3v, cal, colorbar=false, clim = (-1,1), sp =1, title = titles[1],
    color=cgrad(:balance), 
        xlab="σ₁ ≡ m²(Kπ) (GeV²)",
        ylab="σ₃ ≡ m²(pK) (GeV²)")

    cal = [
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
    NaN :
    α2_nonint(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
        
    heatmap!(σ1v, σ3v, cal, colorbar=false, clim = (-1,1), sp = 2, title = titles[2],
    color=cgrad(:balance), 
            xlab="σ₁ ≡ m²(Kπ) (GeV)",
            ylab="σ₃ ≡ m²(pK) (GeV)")

    cal = [
    (Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) > 0 ?
    NaN :
    α3_nonint(Invariants(tbs.ms,σ1=σ1,σ3=σ3); pars=pars)) for σ3 in σ3v, σ1 in σ1v]
                
    heatmap!(σ1v, σ3v, cal, colorbar=true, clim = (-1,1), sp = 3, title = titles[3],
    color=cgrad(:balance), 
            xlab="σ₁ ≡ m²(Kπ) (GeV)",
            ylab="σ₃ ≡ m²(pK) (GeV)")
end
savefig("alpha_heatmap.pdf")



# θhat12 > 0
# θhat23 > 0
# θhat31 > 0
# 
# θhat21 clockwise < 0
# θhat32 clockwise < 0
# θhat13 clockwise < 0

