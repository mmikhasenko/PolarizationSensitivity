
using UpROOT
using Plots
using ThreeBodyDecay
using LinearAlgebra
theme(:wong)

f = TFile(joinpath("C:\\","Users","mikha","HeavyData","Xic2pKpi","2011NtpMagDownMagUp_FinalCuts.root"))
keys(f)
@time t = f["XiR0Tree"][:];
propertynames(t)
length(t)

histogram(t.DalitzX)
histogram2d(t.DalitzX, t.DalitzY, bins=100,
    xlab="m²(Kπ) (GeV²)", ylab="m²(pK) (GeV²)")

