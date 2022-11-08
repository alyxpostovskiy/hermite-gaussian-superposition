begin
    include("main.jl")
    using .HGSuperposition
end

# Define a beam - waist width, wavelength, refractive idnex
beam = Beam(1.0,633e-6, 1.0)

### Definitions for computing superposition coefficients
begin
    # Define an aperture
    sq = Segment([
        [1.0 1.0],
        [-1.0 1.0],
        [-1.0 -1.0],
        [1.0 -1.0]
    ] * 0.5)     # Scale aperture

    # Efield at z=0, power series coefficients
    E0 = [1 0; 0 0;]

    # Highest order mode used
    N = 50

    # Maximum allowable error on E(x,y,0) for any x,y in the aperture
    ϵ = 1e-6
end

# Compute superposition coefficients
C = coef(beam,sq,E0,ϵ,N); 0

# Note: if nothing appears on the first plot, rerun the function.
#   Julia compiles the code on the first run and sometimes nothing gets displayed.
# 2D superposition & plot
begin
    x = -1:0.01:1
    z = 0

    E = superposition(beam,C,x,z); 0
    I = abs2.(E); 0
    plot_2D(I)
end

# 3D superposition & plot
begin
    x = -1:0.02:1
    z = 0:1e2:1e3

    E = superposition(beam,C,x,z); 0
    I = abs2.(E); 0
    plot_3D(x,z,I)

end

# 3D with 2D crossection and intensity cutoff
begin
    x = -1:0.02:1
    z = 0:1e2:1e3

    E = superposition(beam,C,x,z); 0
    I = abs2.(E); 0
    plot_integrated(x,z,I)
end

# We can pass our own coefficient matrix into the superposition
# Complex coefficients work for 2D superpositions
begin
    my_coefs = [0.0 1.0; 
                1.0 0.0]

    x = -2:0.02:2
    z = 0

    E = superposition(beam,my_coefs,x,z); 0
    I = abs2.(E); 0
    plot_2D(I)
end

# Our own coefficients for the 3/2D plot
begin
    x = -2:0.04:2
    z = 0:1e2:1e3

    my_coefs = [0.0 1.0; 
                1.0 0.0]

    E = superposition(beam,my_coefs,x,z); 0
    I = abs2.(E); 0
    plot_integrated(x,z,I,5e-3)
end

# Two sources
begin
    x = -2:0.02:2   # Global picture
    z = 0

    source1_coefs = [1.0 0.0; 0.0 0.0]  #HG_00
    source2_coefs = [-1.0 0.0; 0.0 0.0]  #-HG_00

    # First source has origin (0.5, 0)
    E1 = superposition(beam,source1_coefs,x .- 0.5, x, z); 0
    # Second source has origin (-0.5, 0)
    E2 = superposition(beam,source2_coefs,x .+ 0.5, x, z); 0

    E = E1 + E2; 0
    I = abs2.(E); 0
    plot_2D(I)
end

# Passing in a custom colormap
# Note that to make a colormap / use the predefined ones,
#   you need to explicitly use Makie. This can (and should)
#   be done once at the top of the file with other pkg imports.
using Makie
begin
    my_coefs = [0.0 1.0; 
                1.0 0.0]

    x = -2:0.02:2
    z = 0

    E = superposition(beam,my_coefs,x,z); 0
    I = abs2.(E); 0
    
    # A list of predefined colormaps can be found at
    #   https://docs.makie.org/stable/documentation/colors/
    # You can also define your own from arrays of colors
    plot_2D(I; colormap=:grays)
end