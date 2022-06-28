begin
    include("main.jl")
    using .HGSuperposition
    using CairoMakie
end

# Define an aperture
sq = [
    [1.0 1.0],
    [-1.0 1.0],
    [-1.0 -1.0],
    [1.0 -1.0]
] * 0.5     # Scale aperture

tri = [
    [0.0 1.0],
    [-1.0 -1.0],
    [1.0 -1.0]
] * 0.5

K_reg = K_heuristic(tri,1e-6)

# Calculate coefficients for all modes up to order N
N = 100
K = K_reg(N)
C = coef(N,K,tri); 0

beam = Beam(633e-6, 1.0)
x = -1:0.005:1

E_field = superposition(beam, convert(Array{Float64}, Array(C)), x, 0.0); 0
I = abs2.(E_field); 0

Iplot(I)