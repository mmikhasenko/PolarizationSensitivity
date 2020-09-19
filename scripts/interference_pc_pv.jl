using PolarizationSensitivity
using ThreeBodyDecay
using Statistics
using JLD2
#
using Plots
theme(:wong)
#
const tbs = tbs_Ξc2pKπ
const isobars = (Kst872_pc, Kst872_pv, Δ1232_pc, Δ1232_pv, Λ1520_pc, Λ1520_pv)
const Np = length(isobars)
#

total_pcpv_interference(σs) = sum(real, interference(σs, isobars; i=i, j=j) for i in 1:2:Np, j in 2:2:Np)

plot(tbs.ms, total_pcpv_interference;
    iσx=1, iσy=3, Ngrid=30,
    colorbar=true, c=:balance,
    right_margin=1Plots.PlotMeasures.cm)
