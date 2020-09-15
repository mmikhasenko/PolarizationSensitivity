# 
O(σs,two_λs; pars, isobars) = sum(c*amplitude(σs,two_λs,dc)
    for (c, dc) in zip(pars, isobars))
#
interference(σs, isobars; i,j) = sum(
    conj(O(σs,two_λs; pars = delta(i;N=length(isobars)), isobars=isobars))*
         O(σs,two_λs;pars = delta(j;N=length(isobars)), isobars=isobars) 
    for two_λs in itr(tbs_Ξc2pKπ.two_js))
#
intensity(σs, isobars; pars) = sum(abs2, O(σs,
    ThreeBodySpins(two_h0=two_ν, two_h1=two_λ, two_h2=0, two_h3=0),
    pars=pars, isobars=isobars) for two_λ in [-1,1], two_ν in [-1,1])


parity_map = [1,0,1,1]
function μ_pc(pars; H)
        μ = 0
        for i in 1:length(pars), 
            j in 1:length(pars)
            if parity_map[i]==1 && parity_map[j]==1
                μ += pars'[i]*H[i,j]*pars[j]
            end
        end        
        return μ
end
    
function μ_pv(pars; H)
    μ = 0
    for i in 1:length(pars), 
        j in 1:length(pars)
        if parity_map[i]==0 && parity_map[j]==0
            μ += pars'[i]*H[i,j]*pars[j]
        end
    end        
    return μ
end


