function build_decay_chains(description)
    tbs = description["tbs"]
    chain = description["chain"]
    two_s = description["s"] |> x2
    parity = description["parity"]
    m,Γ = description["mΓ"]
    Ps = description["parities_1230"]
    whichLS = description["whichLS"] # all or min
    # 
    whichLS == "min" && return [decay_chain(chain, (s,σ)->BW(σ, m, Γ);
        two_s = two_s, tbs=tbs, parity = parity, Ps = Ps)]

    whichLS == "all" && return decay_chains(chain, (s,σ)->BW(σ, m, Γ);
        two_s = two_s, tbs = tbs, parity = parity, Ps = Ps)
end

function LSls(description)
    return ThreeBodyDecay.possibleLSls(description["chain"];
        two_js=description["tbs"].two_js,
        two_s = description["s"] |> x2,
        parity = description["parity"],
        Ps=description["parities_1230"])
end

LSls(description)

const tbs = ThreeBodySystem(
    ThreeBodyMasses(m0=2.46867, m1=0.938, m2=0.49367, m3=0.13957),
    ThreeBodySpins(two_h0=1, two_h1=1, two_h2=0, two_h3=0));
    
const tbs_parities_pc = ['+', '-', '-', '+'];
const tbs_parities_pv = ['+', '-', '-', '-'];
# 
build_decay_chains(description)

description = Dict(
    "tbs"=>tbs,
    "chain"=>1,
    "s"=>1,
    "parity"=>'-',
    "mΓ"=>(0.89176, 0.05),
    "parities_1230"=>tbs_parities_pc,
    "whichLS"=>"all")

let σ1 = 0.9^2
    plot()
    for dc in build_decay_chains(description)
        plot!(z->sum(abs2,
            amplitude(
                Invariants(tbs.ms;σ1=σ1,σ3=σ3of1(z,σ1,tbs.ms^2)),
                two_λs,dc)
            for two_λs in itr(tbs.two_js)), -1, 1)
    end
    plot!()
end


fulldescr = [
    Dict(# K*
        "chain"=>1, "s"=>1, "parity"=>'-', "mΓ"=>(0.89176, 0.05),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pc, "whichLS"=>"all"),
    Dict(# Δ
        "chain"=>2, "s"=>3/2, "parity"=>'+', "mΓ"=>(1.232,   0.112),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pc, "whichLS"=>"all"),
    Dict(# Λ
        "chain"=>3, "s"=>3/2, "parity"=>'-', "mΓ"=>(1.5195, 0.017),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pc, "whichLS"=>"all"),
    Dict(# K*
        "chain"=>1, "s"=>1, "parity"=>'-', "mΓ"=>(0.89176, 0.05),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pv, "whichLS"=>"all"),
    Dict(# Δ
        "chain"=>2, "s"=>3/2, "parity"=>'+', "mΓ"=>(1.232,   0.112),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pv, "whichLS"=>"all"),
    Dict(# Λ
        "chain"=>3, "s"=>3/2, "parity"=>'-', "mΓ"=>(1.5195, 0.017),
        "tbs"=>tbs, "parities_1230"=>tbs_parities_pv, "whichLS"=>"all")];
# 

alldcs = vcat(build_decay_chains.(fulldescr)...)[:]

O(σs,two_λs; pars) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(pars, alldcs))

I(σs; pars) = sum(abs2, O(σs,two_λs; pars=pars) for two_λs in itr(tbs.two_js))

sample = flatDalitzPlotSample(tbs.ms; Nev = 10_000);

delta(i; N = 8) = [x==i for x in 1:N];
# 
normpars = [sum(I.(sample; pars=delta(i; N = 8)))/length(sample) for i in 1:8];
normpars .= sqrt.(normpars)
In(σs; pars) = I(σs; pars=pars./normpars)
sum(In.(sample; pars=[0,0,0,1im,0,0,0,0]))/length(sample) ≈ 1.0
#

let
    ps = [plot(σs->In(σs; pars=delta(i)), tbs.ms; iσx = 1, iσy = 3)
        for i in 1:8]
    plot(ps..., layout=grid(2,4), size=(800,300))
end


# 
Interfn(σs,i,j) = sum(
    conj(O(σs,two_λs; pars=delta(i)./normpars)) *
         O(σs,two_λs; pars=delta(j)./normpars)
    for two_λs in itr(tbs.two_js))

H = [sum(Interfn.(sample, i,j)) / length(sample)
        for i in 1:8, j in 1:8]
det(H)

Hext = [[H[i,j] for i in 1:8, j in 1:8] [-1im*H[i,j] for i in 1:8, j in 1:8]
    [H[i,j]*1im for i in 1:8, j in 1:8] [-1im*H[i,j]*1im for i in 1:8, j in 1:8]]
#
rank(Hext)

let B = round.(H; digits=2)
    plot(
        heatmap(real.(H)),
        heatmap(imag.(H)), size=(800,350)
    )
end

plot(σs->real(sum(Interfn(σs,i,j) for i in 1:8, j in 1:8 if i<j)), tbs.ms; iσx = 1, iσy = 3, Ngrid=30)
plot!(colorbar=true)
plot!(color_palette=:balance, clim=(-3,3))

plot(σs->real(Interfn(σs,1,7)), tbs.ms; iσx = 1, iσy = 3)



plot(σs->In(σs; pars=[0,0,0,1,0,0,0,0]), tbs.ms; iσx = 1, iσy = 3)
plot(σs->In(σs; pars=[0,0,0,1,0,0,0,0]), tbs.ms; iσx = 1, iσy = 3)
