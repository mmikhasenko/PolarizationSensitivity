using Plots
using Plots.StatsBase
using PartialWaveFunctions
using ThreeBodyDecay

using TypedTables
using DelimitedFiles
using Cuba
using QuadGK
using Optim
using Test
using ForwardDiff

Np = 4
const tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
#
#                         p    K    π    Ξc
const tbs_parities_pc = ['+', '-', '-', '+'];
const tbs_parities_pv = ['+', '-', '-', '-'];

# chain-1: K*(872)
const Kst872_pc = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Ξc -> K* p
# +  -> -  + # odd wave -> (-1)^L = - => L = 1
# K* J^P = 1^-  ⊗ p J^P = 1/2^+
print(Kst872_pc.two_ls)

const Kst872_pv = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pv)
# Ξc -> K* p
# -  -> -  + # ever wave -> (-1)^L = - => L = 0

# chain-2: Delta**
const Δ1232 = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs,
    parity = '+', Ps = tbs_parities_pc)

# chain-3: Lambda**
const Λ1520  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs,
    parity = '-', Ps = tbs_parities_pc)
# Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
# Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)

delta(i; N=4) = [x==i for x in 1:N]

O(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, (Kst872_pc, Kst872_pv, Λ1520, Δ1232)))

Interf(σs; i,j) = sum(conj(O(σs,two_λs;pars = delta(i))) * O(σs,two_λs;pars = delta(j)) 
    for two_λs in itr(tbs.two_js))

I(σs; pars=[0.1, 0.1, 1.1, 3.9im]) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars) for two_λ in [-1,1], two_ν in [-1,1])
#


#Calculate integrals
Φ = quadgk(σ1->(sqrt(λ(tbs.ms.m0^2, σ1, tbs.ms.m1^2))*sqrt(λ(σ1, tbs.ms.m2^2, tbs.ms.m3^2))/σ1),
     lims1(tbs.ms)...)[1]

const s0 = flatDalitzPlotSample(tbs.ms; Nev = 10_000)

#Calculate all matrix elements Interf_ij
H = Matrix{Complex{Float64}}(undef,Np,Np)

for i in 1:Np
    for j in 1:Np
        H[i,j] = (Φ/length(s0))*sum(Interf.(s0; i=i,j=j))
    end
end

plot(heatmap(log.(abs.(real.(H)))),
     heatmap(imag.(H)))


μ(pars) = real(conj(transpose(pars))*H*pars)

μ(rand(4))

# reading
data = let Nreduced = 100
    M = readdlm(joinpath("data","sims","sample_Kstar=1,1_Lambda=1.1_Delta=3.9.txt"))
    [Invariants(tbs.ms,σ1=M[i,1],σ3=M[i,2]) for i in 1:size(M,1)][1:Nreduced]
end


# let n = 5000   
#     sum(Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms) < 0 
#         for σ1 in 6*rand(n), σ3 in 6*rand(n))/n^2*36
# end

#LIKELIHOOD FUNCTION 
const genpars = [1, 1, 1.1, 3.9]
const genpars′ = genpars./sqrt(μ(genpars)/length(data))
μ(genpars′)
f(genpars′)

ellh(;data,pars) = -sum(log, I(σs; pars) for σs in data) + μ(pars)
#@time ellh(;data = data, pars = rand(4))
#@time I.(data; pars = rand(4))


plot(tbs.ms,σs -> I(σs; pars=[1 + 0im, 1, 1.1, 3.9]),iσx = 1, iσy = 3, c=:viridis,
    colorbar = true, clim = (0,5000))

let Np = 4
    plot()
    for j in 1:4
        ranges = sqrt.(length(data) ./ [H[i,i] for i in 1:4])
    #    init_pars = rand(Np).*ranges
    #    init_pars[1] = rand()*ranges[1]*cis(rand()*2π)
        ph = Vector{Float64}()
        e = Vector{Float64}()
        gen_pars = [1 + 0im, 1, 1.1, 3.9]
        for i in 1:100
            phase = (i-1)/99*2π - π
            modified_pars = copy(gen_pars) 
            modified_pars[j] *= cis(phase)
            append!(ph,phase)
            append!(e,ellh(;data=data,pars = modified_pars))
        end
        plot!(ph,e, lab = "$j")
    #    print(pars)
    #   ellh(;data=data,pars = pars)
    end
    plot!()
end

#Find max values of initial parameters by fulfilling intensity constraint 
#with 1 parameter only
ranges = sqrt.(length(data) ./ [H[i,i] for i in 1:Np])
rand(Np).*ranges.*(cis.(rand(Np)*2π))

fold(x) = x[1:Np] + 1im .* x[(Np+1):2Np]

#Gen init params (in [0,max] range ) with random phase 
init_pars = rand(Np).*ranges.*(cis.(rand(Np)*2π))
print(init_pars[1])
#Scale init params 
init_pars′ = init_pars./sqrt(μ(init_pars)/length(data))
f(x) = ellh(;data=data, pars=x )
f′(x) = fold(ForwardDiff.gradient(p->f(fold(p)), vcat(real(x), imag(x))))
@time f′(init_pars)
f′!(stor,x) = copyto!(stor,f′(x))
Optim.optimize(f, f′!, init_pars′, BFGS(),
                Optim.Options(show_trace = true))


function fit_data!(settings)
    ldata = settings["data"]
    _Natt = settings["Natt"]
    H = settings["H_matrix"]
    # 
    Np = 4
    fold(x) = x[1:Np] + 1im .* x[(Np+1):2Np]
    #
    f(x) = ellh(;data=ldata, pars=x)
    ranges = sqrt.(length(ldata) ./ [H[i,i] for i in 1:Np])
    
    f′(x) = fold(ForwardDiff.gradient(p->f(fold(p)), vcat(real(x), imag(x))))
    f′!(stor,x) = copyto!(stor,f′(x))
    # 
    frs = [let
        init_pars = rand(Np).*ranges.*(cis.(rand(Np)*2π))
        print(init_pars)
        init_pars ./= sqrt(μ(init_pars)/length(data))
        Optim.optimize(f, f′!, init_pars, BFGS(),
                    Optim.Options(show_trace = settings["show_trace"]))
    end for e in 1:_Natt]
    settings["fit_results"] = frs
    print(frs)
end

settings = Dict(
    "H_matrix" => H,
    "data"=> data,
    "Natt"=>3,
    "show_trace"=>false)

fit_data!(settings);


Optim.minimizer(settings["fit_results"][1])

tfr = Table(
    [(st = Optim.converged(fr), min = Optim.minimum(fr), pars = Optim.minimizer(fr))
        for fr in settings["fit_results"]])
print(tfr)

converged = filter(x->x.st, tfr);

print(converged)

histogram(converged.min, bins=100)

bestpars = converged[findmin(converged.min)[2]].pars
print(bestpars)

#const LMs = NamedTuple{(:L,:M)}.([(1,1),(2,1),(4,1)])
#Np = length(LMs)
#fold(x) = x[1:Np] + 1im .* x[(Np+1):2Np]
#init_pars = fold(2rand(2Np))
#print(init_pars)

