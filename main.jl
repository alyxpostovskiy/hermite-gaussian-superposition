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

zr(b::Beam) = pi * b.refr_ind / b.wvlen
k(b::Beam) = 2zr(b)
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

gaussian(a::Complex, x::AbstractFloat) = exp(-a * x^2)

function u_0(q0_q::Complex, x::AbstractFloat)
    norm = √(√(2/pi) * q0_q)
    gauss = gaussian(q0_q, x)
    return norm * gauss
end

function u_1(q0_q::Complex, x::AbstractFloat)
    norm = 2(2/pi)^(1/4)
    gauss = q0_q^(3/2) * x * gaussian(q0_q, x)
    return norm * gauss
end

function u(b::Beam, N::Integer, x::AbstractArray, z::AbstractFloat)
    qz = q(b)(z)
    q0_q = q0(b) / qz
    q̅_q = conj(qz) / qz
    
    xn = length(x)
    ret = Array{ComplexF64}(undef, (N+1, xn))
    # Base case u_0
    for (i, xi) in enumerate(x)
        ret[1, i] = u_0(q0_q, xi)
    end
    # Base case u_1
    if N > 0
        for (i, xi) in enumerate(x)
            ret[2, i] = u_1(q0_q, xi)
        end
    end
    # Recursive case
    if N > 1
        for (i, xi) in enumerate(x), n in 3:(N+1)
            n_dep = 2xi/√(n-1) * q0_q * ret[n-1,i]
            nm1_dep = q̅_q * √((n-2)/(n-1)) * ret[n-2,i]
            ret[n,i] = n_dep + nm1_dep
        end
    end

    return ret
end

function superposition(b::Beam, coefs::AbstractArray, x::AbstractArray, z::AbstractFloat)
    N = size(coefs)[1] - 1
    herms = u(b,N,x,z)
    fourier = exp(-im * k(b) * z)
    field = herms' * coefs * herms
    return field .* fourier
end

"""
Single heatmap plot of an intensity matrix.
"""
function Iplot(I::AbstractMatrix)
    fig = Figure()
    ax = Axis(fig[1,1])
    hidedecorations!(ax)
    heatmap!(ax, I, colormap=:grays)
    return fig
end

end