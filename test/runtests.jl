using Test
using PolarizationSensitivity
using PolarizationSensitivity.ThreeBodyDecay
using QuadGK

@testset "calculation of Dalitz area" begin
    tbs = tbs_Ξc2pKπ

    Φ0 = quadgk(σ1->ρ1(σ1; tbs.ms), lims1(tbs.ms)...)[1]

    # test that integrals are the same
    Φ0_statist = let n = 1000
        sum(Kibble(Invariants(tbs.ms,σ1=σ1,σ3=σ3),tbs.ms^2) < 0 
            for σ1 in 6*rand(n), σ3 in 6*rand(n))/n^2*36
    end
    
    @test abs(Φ0 - Φ0_statist)/Φ0*100 < 10
end


