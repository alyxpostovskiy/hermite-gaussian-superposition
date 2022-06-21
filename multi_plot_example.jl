#=


=#

N = 10:10:100
scale = 0.3:0.1:0.7
x_step = 0.005
x = -1:x_step:1
K = 100

apt = [
    Segment(1.0,-1.0,0.0,2.0)
    Segment(1.0,1.0,-2.0,0.0)
    Segment(-1.0,1.0,0.0,-2.0)
    Segment(-1.0,-1.0,2.0,0.0)
]

N_max = maximum(N)

fig = Figure(resolution=(length(N)*length(x), length(scale)*length(x)), figure_padding=0)
ax = [Axis(fig[i,j]) for i in 1:length(scale), j in 1:length(N)]
hidedecorations!.(ax)
rowgap!(fig.layout,0)
colgap!(fig.layout,0)

a = a_coefs(N_max,K,R_512)

for (i, ρ) in enumerate(scale)
    apt_props = props.(filter(s -> s.dy != 0, apt*ρ))
    print("ψ -- ")
    @time ψ = psi(K, apt_props,R_512)
    print("C -- ")
    @time C = convert(Array{Float64}, Array(coef(a,ψ)))
    print("in -- ")

    A = (2ρ)^2
    @time begin
    for (j, n) in enumerate(N)
        E_field = superposition(beam, C[1:(n+1),1:(n+1)], x, 0.0)
        I = abs2.(E_field)

        xi = Int64(round(200 * (1-ρ)))
        xn = Int64(round(2/x_step))
        I_sums = @view I[(1+xi):(1+xn-xi),(1+xi):(1+xn-xi)]
        μ = sum(I_sums) * x_step^2 / A
        μ2 = sum(I_sums .^ 2) * x_step^2 / A
        σ = sqrt(μ2 - μ^2)

        rowsize!(fig.layout, i, length(x))
        colsize!(fig.layout, j, length(x))
        Label(fig[i,j], "n = $n, ρ = $ρ", color=:white, valign=:top, textsize=36)
        Label(fig[i,j], "μ = $(round(μ,digits=3)), σ = $(round(σ,digits=3))", color=:white, valign=:bottom, textsize=36)
        
        heatmap!(ax[i,j], I, colormap=:grays, colorrange=(0,1.5))
    end
    end
end
fig