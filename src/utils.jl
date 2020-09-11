
delta(i; N) = [x==i for x in 1:N]

arg(z) = atan(imag(z), real(z))

function relative_phase(pars,i,j)
    return arg(pars[i]'*pars[j])
end

function fold(x)
    Np = div(length(x),2)
    return x[1:Np] + 1im .* x[(Np+1):2Np]
end

function contributions(pars; H)  
    return [pars'[i]*H[i,j]*pars[j] for i in 1:length(pars),
        j in 1:length(pars)]  
end

