# tests
# f(x) = ellh(x; data=data, H=H)
# f′ = get_complex_derivative(f)
# myf′!(stor,x) = copyto!(stor,f′(x))

# Optim.optimize(f, myf′!, init_pars, BFGS(),
#                 Optim.Options(show_trace = true))
