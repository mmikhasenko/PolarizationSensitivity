
delta(i; N) = [x==i for x in 1:N]
    
function fold(x)
    Np = div(length(x),2)
    return x[1:Np] + 1im .* x[(Np+1):2Np]
end