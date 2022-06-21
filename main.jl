module HGSuperposition

# Activate project for package versioning
cd(@__DIR__)
using Pkg
Pkg.activate(".")

# Package imports
using Nemo
using CairoMakie
#using PolygonOps

# Export list
export Beam
export K_heuristic
export coef
export superposition
export Iplot

R = ArbField(512)

# Beam type
struct Beam{T<:AbstractFloat}
    wvlen::T
    refr_ind::T
end

k(b::Beam) = 2pi * b.refr_ind / b.wvlen
zr(b::Beam) = k(b) / 2
q0(b::Beam) = zr(b) * im
q(b::Beam) = z -> z + q0(b)
w(b::Beam) = z -> √(1 + (z / zr(b))^2)

# Polygonal Aperture
struct SegmentProps{T<:AbstractFloat}
    υ::T
    ν::T
    y1::T
    y2::T
    function SegmentProps(pt::AbstractArray{S}, dir::AbstractArray{S}) where S <: AbstractFloat
        new{S}(pt[1] - pt[2] * dir[1] / dir[2], dir[1] / (pt[1] * dir[2] - pt[2] * dir[1]), pt[2], pt[2]+dir[2])
    end
end
SegmentProps(t::Tuple) = SegmentProps(t...)

function apt_props(points::AbstractArray)
    ptdir = [(points[i], points[i%length(points) + 1] - points[i]) for i in 1:length(points)]
    nonhoriz = filter(x -> x[2][2] != zero(x[2][2]), ptdir)
    return SegmentProps.(nonhoriz)
end

# Vec to Diag Mat
diagonal(v,n) = [i==j ? v[i] : 0 for i in 1:n, j in 1:n]

# Taylor series coefs
function a_coefs(N::Integer,K::Integer,R::Nemo.ArbField)
    if N%2 == 1 || K%2 == 1
        error("N and K must be even")
    end

    a = Nemo.zero_matrix(R,N+1,K+1)
    a_coefs!(a)
    return a
end

function a_coefs!(zero_mat::arb_mat)
    R = base_ring(zero_mat)
    N, K = size(zero_mat)
    N -= 1; K -= 1

    norm = √(√(R(2)/R(pi)))

    for k in 0:(K÷2)    #a_{0,k}
        zero_mat[1,2k+1] = norm * (-1)^k / R(factorial(big(k)))
    end
    for k in 1:(K÷2)    #a_{1,k}
        zero_mat[2,2k] = -2k * zero_mat[1,2k+1]
    end
    for n in 1:(N÷2)    #a_{n,0}
        zero_mat[2n+1,1] = (-1)^n * norm * √(R(factorial(big(2n)))) / R(big(2)^n * factorial(big(n)))
    end
    for k in 1:(K÷2), n in 1:(N-1)
        zero_mat[n+2,2k+n%2] = 2/√(R(n+1)) * zero_mat[n+1,2k-1+n%2] - √(R(n)/R(n+1)) * zero_mat[n,2k+n%2]
    end
    return zero_mat
end

# Aperture characteristic
function psi(K::Integer,apt::AbstractVector,R::Nemo.ArbField)
    ψ = Nemo.zero_matrix(R,K+1,K+1)
    seg_props = apt_props(apt)
    S = length(seg_props)


    factorials = [R(factorial(i)) for i in big(0):big(K+1)]
    rec_facs = R(1) ./ factorials
    nCr = Nemo.matrix(R, [j > i ? R(0) : factorials[i] * rec_facs[i-j+1] * rec_facs[j+1] for i in 1:(K+1), j in 0:(K+1)])

    υ = [(R(s.υ)^i) for s in seg_props, i in 1:(K+1)]
    ν = [R(s.ν)^i for s in seg_props, i in 0:(K+1)]
    y = [(R(s.y2)^i - R(s.y1)^i)/i for s in seg_props, i in 1:(2K+2)]

    for s in 1:S
        υs = Nemo.matrix(R, diagonal(υ[s,:],K+1))
        νs = Nemo.matrix(R, diagonal(ν[s,:],K+2))
        Γ = υs * nCr * νs
        Y = Nemo.matrix(R, [y[s,k+j] for k in 0:(K+1), j in 1:(K+1)])

        ψ += Γ * Y
    end

    return ψ
end

# What K is needed for all n in N to be calculated to relative error ϵ
function K_terms_rel(a::Nemo.arb_mat,ψ::Nemo.arb_mat,N::AbstractVector,ϵ::AbstractFloat) # a is r x k, ψ is k x k, all n, k even
    ret = zeros(Int64, length(N))
    i = 1
    C_old = (a[(N[i]+1):(N[i]+1),1:1]) * (ψ[1:1,1:1]) * Nemo.transpose(a[(N[i]+1):(N[i]+1),1:1])
    for k in 2:(size(a)[2]-1)÷2
        C_cur = (a[(N[i]+1):(N[i]+1),1:2k+1]) * (ψ[1:2k+1,1:2k+1]) * Nemo.transpose(a[(N[i]+1):(N[i]+1),1:2k+1])
        err = 2 * abs((C_cur[end,end] - C_old[end,end]) / (C_cur[end,end] + C_old[end,end]))
        if err < ϵ
            ret[i] = 2k
            if i == length(N)
                return ret
            end
            i += 1
        end
        C_old = C_cur
    end
    return ret
end

"""
Given an aperture, a maximum error, and (optionally) an initial K,
returns a regressed function of n, which provides a K value for that
N for which the error in the coefficient is less than ϵ.
"""
function K_heuristic(apt,ϵ::AbstractFloat,K_0::Integer=150)
    N_0 = 500
    R_512 = ArbField(512)
    a = a_coefs(N_0,K_0,R_512)
    ψ = psi(K_0,apt,R_512)

    N = 100:100:500
    nv = [ones(length(N)) (N .^ (2/3))]
    kv = K_terms_rel(a,ψ,N,ϵ)
    if any(map(x -> x == zero(x), kv))
        error("Increase K_0 and try again")
    end

    reg = nv\kv

    K_reg(x) = 2*Int64(ceil(0.5*(reg[1] + reg[2] * (x^(2/3)))))
    return K_reg
end

coef(a::arb_mat,ψ::arb_mat) = a * ψ * transpose(a)
coef(a::arb_mat,ψ::arb_mat,N::Integer,K::Integer) = coef(a[1:(N+1),1:(K+1)],ψ[1:(K+1),1:(K+1)])
function coef(N::Integer,K::Integer,apt::AbstractVector)
    a = a_coefs(N,K,R)
    ψ = psi(K,apt,R)
    return coef(a,ψ)
end

# calculates herm polys for all degrees up to n
function h_J(arg::AbstractFloat, phase::Number, n::Integer)
    h = Array{typeof(phase)}(undef, n + 1)
    h[1] = one(phase)
    if n > 0
        h[2] = phase * arg
    end
    if n > 1
        for i in 2:n
            h[i+1] = (phase / √i) * arg * h[i] - √((i - 1) / i) * h[i-1]
        end
    end
    return h
end

function superposition(b::Beam, coefs::AbstractArray, x::AbstractArray, z::AbstractFloat)
    n = size(coefs)[1] - 1

    qb = q(b)
    wb = w(b)
    kb = k(b)

    xn = length(x)
    
    phase = √(-conj(qb(z)) / qb(z))
    herm_arg = @. (2 / wb(z)) * x
    herms = Array{typeof(phase)}(undef, n+1, xn)
    for xi in 1:xn
        herms[:,xi] = h_J(herm_arg[xi], phase, n)
    end

    field = herms' * coefs * herms

    gauss = @. exp((-im * kb) / (2*qb(z)) * x^2)
    gauss_xy = reshape(gauss, (xn,1)) .* reshape(gauss, (1,xn))
    field .*= gauss_xy

    z_factor = √(2/pi) * (q0(b)/qb(z)) * exp(-im * k(b) * z)

    return field .* z_factor
end

"""
Single heatmap plot of an intensity matrix.
"""
function Iplot(I::AbstractMatrix)
    fig = Figure()
    ax = Axis(fig[1,1])
    hidedecorations!(ax)
    heatmap!(ax, I, colormap=:grays, colorrange=(0,1.5))
    return fig
end

end