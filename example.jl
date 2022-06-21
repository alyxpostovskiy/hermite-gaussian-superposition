include("main.jl")
using .HGSuperposition
using CairoMakie

# Define an aperture
sq = [
    [1.0 1.0],
    [-1.0 1.0],
    [-1.0 -1.0],
    [1.0 -1.0]
] * 0.5     # Scale aperture

K_reg = K_heuristic(sq,1e-6)

# Calculate coefficients for all modes up to order N
N = 1000
K = K_reg(N)
C = coef(N,K,sq)

beam = Beam(633e-6, 1.0)
x = -1:0.005:1

E_field = superposition(beam, convert(Array{Float64}, Array(C)), x, 0.0)
I = abs2.(E_field)

Iplot(I)