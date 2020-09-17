using Test
using PolarizationSensitivity
using PolarizationSensitivity.ThreeBodyDecay
using QuadGK


tbs = tbs_Ξc2pKπ

@testset "calculation of Dalitz area" begin

    Φ0 = quadgk(σ1->ρ1(σ1; tbs.ms), lims1(tbs.ms)...)[1]

    # test that integrals are the same
    Φ0_statist = let n = 5000
        sum(Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) < 0 
            for σ1 in 6*rand(n), σ3 in 6*rand(n))/n^2*36
    end
    
    @test abs(Φ0 - Φ0_statist)/Φ0*100 < 10
end


@testset "#Test that αs return a value between -1,1" begin
    isobars = (Kst872_pc, Kst872_pv)
    sample = flatDalitzPlotSample(tbs.ms; Nev = 10)
    pars = [1,1]
    α1 =  α1_avg(pars; isobars = isobars, sample = sample)
    @test -1 <= α1 <= 1
    α2 =  α2_avg(pars; isobars = isobars, sample = sample)
    @test -1 <= α2 <= 1
    α3 =  α3_avg(pars; isobars = isobars, sample = sample)
    @test -1 <= α3 <= 1 
end

