

const tbs_parities_pc = ['+', '-', '-', '+'];
const tbs_parities_pv = ['+', '-', '-', '-'];

# chain-1: K*(872)
# Ξc -> K* p
# +  -> -  + # odd wave -> (-1)^L = - => L = 1
# K* J^P = 1^-  ⊗ p J^P = 1/2^+
const Kst872_pc = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs_Ξc2pKπ,
    parity = '-', Ps = tbs_parities_pc)

# Ξc -> K* p
# -  -> -  + # ever wave -> (-1)^L = - => L = 0
const Kst872_pv = decay_chain(1, (s,σ)->BW(σ, 0.89176, 0.05);
    two_s = 1|>x2, tbs=tbs_Ξc2pKπ,
    parity = '-', Ps = tbs_parities_pv)

# chain-2: Delta**
const Δ1232_pc = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs_Ξc2pKπ,
    parity = '+', Ps = tbs_parities_pc)
# 
const Δ1232_pv = decay_chain(2, (s,σ)->BW(σ, 1.232,   0.112);
    two_s = 3/2|>x2, tbs=tbs_Ξc2pKπ,
    parity = '+', Ps = tbs_parities_pv)

# chain-3: Lambda**
const Λ1520_pc  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs_Ξc2pKπ,
    parity = '-', Ps = tbs_parities_pc)
#

const Λ1520_pv  = decay_chain(3, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs_Ξc2pKπ,
    parity = '-', Ps = tbs_parities_pv)
# Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
# Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)
