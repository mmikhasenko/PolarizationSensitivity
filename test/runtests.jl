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


@testset "αs return a value between -1,1" begin
    isobars = (Kst872_pc, Kst872_pv)
    sample = flatDalitzPlotSample(tbs.ms; Nev = 10)
    pars = [1,1]
    α1 =  α_avg(I1_assym, pars; isobars = isobars, sample = sample)
    @test -1 <= α1 <= 1
    α2 =  α_avg(I2_assym,pars; isobars = isobars, sample = sample)
    @test -1 <= α2 <= 1
    α3 =  α_avg(I3_assym,pars; isobars = isobars, sample = sample)
    @test -1 <= α3 <= 1 
end

@testset "parsing datafile name" begin
    datafile = joinpath("data","sims","sample_Kstar=1.3,1-1im_Delta=2-0.6im,2+1im_Lambda=1.2-0.5im,2+0.3im.txt")
    # 
    @test_throws ErrorException parse_values_from_datafile_name(datafile, 3) # Np =! 6
    # 
    genpars = parse_values_from_datafile_name(datafile, 6)
    @test genpars ≈ [1.3, 1-1im, 2-0.6im, 2+1im, 1.2-0.5im,2+0.3im]
end

