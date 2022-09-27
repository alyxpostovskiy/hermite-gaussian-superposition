module HGSuperposition

# Activate project for package versioning
begin
    cd(@__DIR__)
    using Pkg
    Pkg.activate(".")
end

# Package imports
using Nemo
using WGLMakie
using JSServe
using Observables

# Exports
begin
    export Beam
    export Segment
    export coef
    export superposition
    export plot_2D
    export plot_3D
    export plot_integrated
end

# Default globals
const R_def = Nemo.ArbField(512)
const cmap_def = :thermal

# Beam type
begin
    """Gaussian beam with parameters waist width, 
        wavelength, and index of refraction of medium."""
    struct Beam{T<:AbstractFloat}
        w_0::T
        wvlen::T
        refr_ind::T
    end
    zr(b::Beam) = pi * b.refr_ind / b.wvlen
    k(b::Beam) = 2zr(b)
    q0(b::Beam) = zr(b) * im
    q(b::Beam) = z -> z + q0(b)
    w(b::Beam) = z -> √(1 + (z / zr(b))^2)
end

# Aperture
begin
    ## An aperture is a vector of row matrices [x y]
    """A line segment from pt1 to pt2 in R^2.
        Can be constructed from two points or a tuple of two points.
        Passing in a vector representing a continuous path will return
        a vector of the path segments; no check to ensure the path is valid is performed."""
    struct Segment{T<:AbstractFloat}
        x1::T
        y1::T
        x2::T
        y2::T
        function Segment(pt1::AbstractArray{S}, pt2::AbstractArray{S}) where S <: AbstractFloat
            new{S}(pt1[1], pt2[1], pt1[2], pt2[2])
        end
    end
    Segment(t::Tuple) = Segment(t...)
    function Segment(points::AbstractVector)
        pt_pairs = [(points[i], points[i%length(points) + 1]) for i in eachindex(points)]
        return Segment.(pt_pairs)
    end

    nonhoriz_segs(segs::AbstractArray{T}) where T <: Segment = filter(s -> s.y1 != s.y2, segs)

    seg_det(s::Segment) = s.x1 * s.y2 - s.y1 * s.x2
    seg_υ(s::Segment) = seg_det(s)/(s.y2 - s.y1)
    seg_ν(s::Segment) = (s.x2 - s.x1) / seg_det(s)

    # Furthest point
    apt_max(segs::AbstractArray{T}) where T <: Segment = maximum([max(abs(s.x1),abs(s.y1)) for s in segs])

    # By shoelace formula
    apt_area(segs::AbstractArray{T}) where T <: Segment = sum(seg_det.(segs))/2
end

# Matrix generation and util
begin
    # Vec of length N2-N1+1 where v[n] = x^(N1+n-1)
    pow_vec(x,N::Integer) = pow_vec(x,0,N)
    pow_vec(x,N1::Integer,N2::Integer) = [x^n for n in N1:N2]

    # Vec to Diag Mat
    diagonal(v) = [i==j ? v[i] : 0 for i in eachindex(v), j in eachindex(v)]

    # Diagonal N+1 x N+1 matrix where D[n,n] = x^(n-1)
    pow_dmat(x,N::Integer) = diagonal(pow_vec(x,N))
    pow_dmat(x,N1::Integer,N2::Integer) = diagonal(pow_vec(x,N1,N2))
    pow_dmat(x,N::Integer,R::Nemo.ArbField) = Nemo.matrix(R, pow_dmat(x,N))
    pow_dmat(x,N1::Integer,N2::Integer,R::Nemo.ArbField) = Nemo.matrix(R, pow_dmat(x,N1,N2))

    pow_int(x1,x2,p) = (x2^(p+1) - x1^(p+1))/(p+1)
    pow_int_mat(x1,x2,N::Integer,M::Integer) = [pow_int(x1,x2,i+j-2) for i in 1:(N+1), j in 1:(M+1)]
    pow_int_mat(x1,x2,N::Integer,M::Integer,R::Nemo.ArbField) = Nemo.matrix(R, pow_int_mat(x1,x2,N,M))

    #N+1 by N+2 matrix
    Γ_gen(N::Integer,R::Nemo.ArbField) = 
        Nemo.matrix(R, [i > j - 1 ? 
                        factorial(R(i-1)) / (factorial(R(j-1)) * factorial(R(i - j + 1))) : 
                        0 for i in 1:(N+1), j in 1:(N+2)])
end


# Computing aperture characteristic
begin
    function T_gen(b::Beam, apt::AbstractVector{T}, N::Integer, M::Integer, R::Nemo.ArbField=R_def) where T <: Segment
        Γ_N = Γ_gen(N,R)
        return sum(pow_dmat(seg_υ(s)/b.w_0,1,N+1,R) * Γ_N * pow_dmat(seg_ν(s)*b.w_0,N+1,R) * pow_int_mat(s.y1,s.y2,N+1,M,R) for s in nonhoriz_segs(apt))
    end

    function T_gen(b::Beam, apt::AbstractVector{T}, K::Integer, E0::AbstractMatrix, R::Nemo.ArbField=R_def) where T <: Segment
        return T_gen(b,apt,K+size(E0)[1]-1,K+size(E0)[2]-1,R)
    end

    function P_gen(beam::Beam, T, E0)
        K = min(size(T)[1] - size(E0)[1], size(T)[2] - size(E0)[2])
        return sum(E0[i,j] * beam.w_0^(i+j-2) * T[i:i+K,j:j+K] for i in 1:size(E0)[1], j in 1:size(E0)[2])
    end

    function P_gen(beam::Beam, apt::AbstractVector{T}, K::Integer, E0::AbstractMatrix, R::Nemo.ArbField=R_def) where T <: Segment
        return P_gen(beam, T_gen(beam,apt,K,E0,R), E0)
    end
end


# Taylor series coefs and error bound
begin
    function a_nk_direct(n::Integer,k::Integer)
        norm = √(√(2/π) * factorial(n) / 2^n) * (-1)^((n-k)/2)
        return norm * sum((2√2)^(k-2i) / (factorial(k-2i) * factorial(Integer((n-k)/2) + i) * factorial(i)) for i in max(0,Integer((k-n)/2)):Integer(floor(k/2)))
    end

    #=
    function K_terms(ϵ,L,n::Integer)
        if ϵ < 0
            throw(DomainError(ϵ, "ϵ must be greater than 0"))
        end

        ϵ_trans = π^(-1/4) - √(π^(-1/2) - ϵ/4)
        k = n%2
        while true
            if abs(a_nk_direct(big(n),big(k+2))) * abs(L)^(k+3) / (k+3) < ϵ_trans
                return k
            end
            k+=2
        end
    end
    =#

    function K_terms(ϵ,L,A,E0::AbstractMatrix,n::Integer)
        if ϵ < 0
            throw(DomainError(ϵ, "ϵ must be greater than 0"))
        end

        ϵ_trans = ϵ / (4A * π^(-1/2) * (n+1)^2)
        k = n%2
        while true
            bound = abs(a_nk_direct(big(n),big(k+2))) * 
                    sum(abs(E0[i,j]) * abs(L)^(i+j+k+2) *
                        (1/((i+k+2)*(j)) + 1/((i)*(j+k+2)))
                        for i in 1:size(E0)[1], j in 1:size(E0)[2])
            if bound < ϵ_trans
                return k
            end
            k+=2
        end
    end

    function K_terms(ϵ::Real,apt::AbstractArray{T},E0::AbstractMatrix,n::Integer) where T<:Segment
        L = apt_max(apt)
        A = apt_area(apt)
        return K_terms(ϵ,L,A,E0,n)
    end

    # Taylor series coefs
    function a_coefs(N::Integer,K::Integer,R::ArbField=R_def)
        a = Nemo.zero_matrix(R,N+1+N%2,K+1+K%2)
        a_coefs!(a)
        return a[1:end-N%2,1:end-K%2]
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
end


# Coefficient computation
begin
    """Computes a matrix of superposition coefficients for a given
        beam, aperture, electric field at the waist, maximum allowable error for E(x,y,0),
        and the highest order of mode to compute the coefficients for 
        (i.e. N=4 will compute 25 coefficients, c_00 through c_44).
        Returns an arb_matrix."""
    coef(A,P) = A * P * transpose(A)
    function coef(b::Beam,apt::AbstractVector{T}, E0::AbstractMatrix, ϵ::Real, N::Integer, R::ArbField=R_def) where T <: Segment
        K = K_terms(ϵ,apt,E0,N)
        A = a_coefs(N,K,R)
        P = P_gen(b,apt,K,E0,R)
        return coef(A,P)
    end
end


# Superposition
begin
    gaussian(a::Complex, x::Real) = exp(-a * x^2)

    function u_0(q0_q::Complex, x::Real)
        norm = √(√(2/pi) * q0_q)
        gauss = gaussian(q0_q, x)
        return norm * gauss
    end

    function u_1(q0_q::Complex, x::Real)
        norm = 2(2/pi)^(1/4)
        gauss = q0_q^(3/2) * x * gaussian(q0_q, x)
        return norm * gauss
    end

    function u(b::Beam, N::Integer, x::AbstractArray, z::Real)
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

    """Computes a superposition of Hermite-Gaussian modes for a given
    Gaussian beam, coefficient matrix (real or complex) over the given coordinates.
    Values of the superposition will be computed on the Cartesian product of x, x, and z.
    Note that z can be a number, in which case a 2D superposition will be returned."""
    function superposition(b::Beam, coefs::AbstractArray, x::AbstractArray, z::Real)
        N = size(coefs)[1] - 1
        herms = u(b,N,x,z)
        fourier = exp(-im * k(b) * z)
        field = herms' * coefs * herms
        return field .* fourier
    end

    function superposition(b::Beam, coefs::AbstractArray, x::AbstractArray, z::AbstractArray)
        return permutedims(cat((superposition(b,coefs,x,zi) for (i,zi) in enumerate(z))...; dims=3), [3,1,2])
    end

    superposition(b::Beam,coefs::arb_mat,x,z) = superposition(b,convert(Array{Float64},Array(coefs)),x,z)
end


# Plotting Util
begin
    slider_DOM(slider, text) = DOM.div(text, slider, slider.value)
    function slider_DOM(slider, text, val_map)
        vmap_app = map(slider) do idx
            return val_map(idx)
        end
        return DOM.div(text, slider, vmap_app)
    end

    function plot_DOM(dom)
        return App() do session::Session
            return JSServe.record_states(session, dom)
        end
    end
end

# Plotting
begin
    """Creates a 2D heatmap plot of an intensity matrix.
        A colormap can be provided; by default, :thermal is used."""
    function plot_2D(I;colormap=cmap_def)
        fig = Figure()
        ax = Axis(fig[1,1])
        hidedecorations!(ax)
        heatmap!(ax, I, colormap=colormap)
        return fig
    end

    """Creates a 3D Maximum Intensity Projection of an intensity array.
        A scale for z can be provided; by default the plot will be scaled such
            that the z-axis is as long as one unit on the xy-axes.
        A colormap can be provided; by default, :thermal is used."""
    function plot_3D(x,z,I, z_scale::Real; colormap=cmap_def)
        fig = Figure()
        lscene = LScene(fig[1,1], show_axis=false, scenekw=(clear=true,backgroundcolor=:white,))
        I_norm = I ./ maximum(I)

        volume!(lscene, z .* z_scale, x, x, I_norm, algorithm=:mip, colormap=colormap)

        return plot_DOM(DOM.div(fig))
    end
    plot_3D(x,z,I;colormap=cmap_def) = plot_3D(x,z,I,1/maximum(z);colormap)

    """Creates a 3D Maximum Intensity Projection of an intensity array with
        an accompanying 2D heatmap of a cross-section. The z-position of the
        cross-section varies with a slider. An additional slider functions as
        an intensity cutoff, such that all intensities above the slider value are
        considered equal by the plot.
        A scale for z can be provided; by default the plot will be scaled such
            that the z-axis is as long as one unit on the xy-axes.
        A colormap can be provided; by default, :thermal is used.
        A resolution for the whole plot (in pixels) can be provided;
            by default, it is (1200,800)."""
    function plot_integrated(x,z,I,z_scale::Real; colormap=cmap_def, resolution=(1200,800))
        fig = Figure(resolution=resolution)
        lscene = LScene(fig[1,1], show_axis=false, scenekw=(clear=true,backgroundcolor=:white,))
        ax = Axis(fig[1,2])
        hidedecorations!(ax)

        sliders = Dict("z" => Slider(1:length(z)), "cutoff" => Slider(0:0.05:1))
        sliders_DOM = [slider_DOM(sliders["z"], "z: ", i -> z[i]), 
                        slider_DOM(sliders["cutoff"], "cutoff: ")]

        I_norm = I ./ maximum(I)
        slice = map(sliders["z"]) do idx
            return I_norm[idx,:,:]
        end
        I_plot = map(sliders["cutoff"]) do idx
            return I_norm ./ idx
        end

        volume!(lscene, z .* z_scale, x, x, I_plot, algorithm=:mip, colormap=colormap)
        heatmap!(ax, slice, colormap=:thermal, colorrange=(0,1))

        return plot_DOM(DOM.div(fig, sliders_DOM; style="color: red"))
    end
    plot_integrated(x,z,I;colormap=cmap_def, resolution=(1200,800)) = plot_integrated(x,z,I,1/maximum(z);colormap,resolution)
end

end