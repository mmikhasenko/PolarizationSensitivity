
μ(pars; H) = real(pars' * H * pars)

ellh(pars,isobars;data,H) = -sum(log, intensity(σs, isobars; pars=pars) for σs in data) + μ(pars;H=H)

#
function get_complex_derivative(f)
    f′(x) = fold(ForwardDiff.gradient(p->f(fold(p)), vcat(real(x), imag(x))))
    return f′
end

function get_inplace_derivative(f)
    f′ = get_complex_derivative(f)
    f′!(stor,x) = copyto!(stor,f′(x))
    return f′!
end

function fit_data!(settings)
    ldata = settings["data"]
    _Natt = settings["Natt"]
    H = settings["H_matrix"]
    # 
    f(x) = ellh(x;data=ldata, H=H)
    f′ = get_complex_derivative(f)
    f′!(stor,x) = copyto!(stor, f′(x))
    # 
    # Find max values of initial parameters by fulfilling intensity constraint 
    # with 1 parameter only
    ranges = sqrt.(length(ldata) ./ [H[i,i] for i in 1:Np])
    #
    frs = [let
        # Gen init params (in [0,max] range ) with random phase 
        init_pars = rand(Np).*ranges.*(cis.(rand(Np)*2π))
        init_pars ./= sqrt(μ(init_pars; H=H)/length(data))
        Optim.optimize(f, f′!, init_pars, BFGS(),
                    Optim.Options(show_trace = settings["show_trace"]))
    end for e in 1:_Natt]
    settings["fit_results"] = frs
end
