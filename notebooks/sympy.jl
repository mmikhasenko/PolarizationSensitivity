### A Pluto.jl notebook ###
# v0.11.2

using Markdown
using InteractiveUtils

# ╔═╡ 76a8b910-d717-11ea-1f02-4d6a1f325869
md"""
## Dependences
Call the symbolic Wigner D function and Clebsch-Gordan coefficient from the `sympy.physics` module.
"""

# ╔═╡ 2a032410-d717-11ea-020d-33211d6293b1
md"""
### Convenient signature of the functions
"""

# ╔═╡ b282fd10-d717-11ea-2199-6d3e5bb43e82
md"""
### Start coding and symplifying expressions
"""

# ╔═╡ 3b809190-d718-11ea-21c7-7b8626476448
md"""
The following expression is simplified
``
$
V_\lambda(\theta) = \sum_{\xi}^{[-1,1]} \left( d_{\lambda,\xi}^{1}(\theta) \right)^2
$
``
"""

# ╔═╡ 05882a00-d71b-11ea-2c5f-bde3a3b579b2
md"""
The expressions are transformed to an evaluatable code and plotted below.
"""

# ╔═╡ a6f0c230-d716-11ea-0545-f7a12b7334d0
begin
	using SymPy
	#
	import PyCall
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:WignerD, :wignerd), typ=:Any)
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	import_from(sympy.physics.wigner)
end

# ╔═╡ 63457f20-d717-11ea-20af-e7566bebb3b4
using LinearAlgebra

# ╔═╡ 70bb6390-d717-11ea-3e41-af29b6d47470
using Parameters

# ╔═╡ f15aa290-d717-11ea-3835-8f7428eb9d20
using Plots

# ╔═╡ fbc00820-d716-11ea-1777-1d0a08d57a83
Wignerd(J,λ1,λ2,θ) = WignerD(J,λ1,λ2,0,θ,0).doit()

# ╔═╡ ace19a20-d716-11ea-2eee-83d1db602f9e
clgn(two_j1,two_m1,two_j2,two_m2,two_j,two_m) = clebsch_gordan(
		Sym(two_j1)/2, Sym(two_j2)/2, Sym(two_j)/2,
		Sym(two_m1)/2, Sym(two_m2)/2, Sym(two_m)/2)

# ╔═╡ b109e710-d716-11ea-1a6b-fd3c9dc53a16
θ,ϕ = @vars θ ϕ real=true

# ╔═╡ b5caca80-d716-11ea-31d8-1149822f3ef9
begin
	V(λ) = sum(abs2, Wignerd(1,λ,ξ,θ) for ξ in [-1, 1])
	v = simplify.([V(λ) for λ = [-1, 0, 1]])
end

# ╔═╡ d9201700-d717-11ea-19d2-19e2746dd7ed
expr = lambdify.(v)

# ╔═╡ e97a2fa0-d717-11ea-29ce-b1dac4f950a2
let
	plot(ylab="Vλ(θ)", xlab="θ")
	for (i,e) in enumerate(expr)
		plot!(e, -π/2, π, lab="λ=$(i-2)")
	end
	plot!()
end

# ╔═╡ Cell order:
# ╟─76a8b910-d717-11ea-1f02-4d6a1f325869
# ╟─a6f0c230-d716-11ea-0545-f7a12b7334d0
# ╟─2a032410-d717-11ea-020d-33211d6293b1
# ╠═fbc00820-d716-11ea-1777-1d0a08d57a83
# ╠═ace19a20-d716-11ea-2eee-83d1db602f9e
# ╟─b282fd10-d717-11ea-2199-6d3e5bb43e82
# ╠═b109e710-d716-11ea-1a6b-fd3c9dc53a16
# ╟─3b809190-d718-11ea-21c7-7b8626476448
# ╠═b5caca80-d716-11ea-31d8-1149822f3ef9
# ╟─05882a00-d71b-11ea-2c5f-bde3a3b579b2
# ╠═d9201700-d717-11ea-19d2-19e2746dd7ed
# ╠═e97a2fa0-d717-11ea-29ce-b1dac4f950a2
# ╠═63457f20-d717-11ea-20af-e7566bebb3b4
# ╠═70bb6390-d717-11ea-3e41-af29b6d47470
# ╠═f15aa290-d717-11ea-3835-8f7428eb9d20
