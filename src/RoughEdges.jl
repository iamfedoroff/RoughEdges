module RoughEdges

import FFTW
import Images
import Interpolations: linear_interpolation
import KernelAbstractions: @index, @kernel, get_backend
import Random

export rough, geometry_mask


# ******************************************************************************************
# Power spectral density (PSD) functions
#
# C.A. Mack, "Analytic form for the power spectral density in one, two, and three
# dimensions", Journal of Micro/Nanolitography, MEMS, and MOEMS, 10, 040501 (2011)
# ******************************************************************************************
function psd_exp(f, sigma, xi)
    return 2 * sigma^2 * xi / (1 + (2*pi * f * xi)^2)
end

function psd_exp(fx, fy, sigma, xix, xiy)
    return 2*pi * sigma^2 * xix * xiy /
           (1 + (2*pi * fx * xix)^2 + (2*pi * fy * xiy)^2)^(3/2)
end

function psd_exp(fx, fy, fz, sigma, xix, xiy, xiz)
    return 8*pi * sigma^2 * xix * xiy * xiz /
           (1 + (2*pi * fx * xix)^2 + (2*pi * fy * xiy)^2 + (2*pi * fz * xiz)^2)^2
end


function psd_gauss(f, sigma, xi)
    return sqrt(pi) * sigma^2 * xi * exp(-(pi * f * xi)^2)
end

function psd_gauss(fx, fy, sigma, xix, xiy)
    return pi * sigma^2 * xix * xiy * exp(-(pi * fx * xix)^2) * exp(-(pi * fy * xiy)^2)
end

function psd_gauss(fx, fy, fz, sigma, xix, xiy, xiz)
    return pi^(3/2) * sigma^2 * xix * xiy * xiz *
           exp(-(pi * fx * xix)^2) * exp(-(pi * fy * xiy)^2) * exp(-(pi * fz * xiz)^2)
end


# ******************************************************************************************
# Rough edge, surface, and volume generators
#
# C.A. Mack, "Generating random rough edges, surfaces, and volumes", Applied Optics, 52,
# 1472 (2013); https://doi.org/10.1364/AO.52.001472
#
# Spectral grid:
#       f[1] - highest frequency
#   f[2:N/2] - negative frequencies
#   f[N/2+1] - zero frequency
# f[N/2+2:N] - positive frequencies
# ******************************************************************************************
function rough(xin; ac=:gauss, sigma=nothing, xi=nothing, seed=nothing)
    Nxin = length(xin)
    isodd(Nxin) ? x = evenize(xin) : x = xin

    Nx = length(x)
    dx = x[2]-x[1]
    Lx = x[end]-x[1]

    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))

    if ac in (:exp, :gauss)
        if isnothing(sigma) || isnothing(xi)
            error("For '$ac' autocorrelation you have to specify 'sigma' and 'xi'")
        end
        if ac == :exp
            psd = psd_exp
        elseif ac == :gauss
            psd = psd_gauss
        end
        PSD = [psd(fxn, sigma, xi) for fxn=fx]
    elseif typeof(ac) <: Function
        PSD = [ac(xn) for xn=x]
        PSD = FFTW.ifftshift(FFTW.fft(FFTW.fftshift(PSD))) * dx
    end

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx)
    for ix=1:Nx
        if ix == 1   # highest frequency
            F[ix] = sqrt(Lx * PSD[ix]) * randn()
        elseif fx[ix] == 0   # zero frequency
            F[ix] = sqrt(Lx * PSD[ix]) * randn()
        elseif fx[ix] > 0   # positive frequencies
            F[ix] = sqrt(Lx * PSD[ix]) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for ix=2:Nx
        if fx[ix] < 0   # negative frequencies
            F[ix] = conj(F[Nx - ix + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / dx

    isodd(Nxin) ? R = R[1:end-1] : nothing

    check_imag(R)

    return real.(R)
end


function rough(xin, yin; ac=:gauss, sigma=nothing, xi=nothing, seed=nothing)
    Nxin, Nyin = length(xin), length(yin)
    isodd(Nxin) ? x = evenize(xin) : x = xin
    isodd(Nyin) ? y = evenize(yin) : y = yin

    Nx, Ny = length(x), length(y)
    dx, dy = x[2]-x[1], y[2]-y[1]
    Lx, Ly = x[end]-x[1], y[end]-y[1]

    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))
    fy = FFTW.ifftshift(FFTW.fftfreq(Ny, 1/dy))

    if ac in (:exp, :gauss)
        if isnothing(sigma) || isnothing(xi)
            error("For '$ac' autocorrelation you have to specify 'sigma' and 'xi'")
        end
        if ac == :exp
            psd = psd_exp
        elseif ac == :gauss
            psd = psd_gauss
        end
        if typeof(xi) <: Real
            xix = xiy = xi
        else
            xix, xiy = xi
        end
        PSD = [psd(fxn, fyn, sigma, xix, xiy) for fxn=fx, fyn=fy]
    elseif typeof(ac) <: Function
        PSD = [ac(xn, yn) for xn=x, yn=y]
        PSD = FFTW.ifftshift(FFTW.fft(FFTW.fftshift(PSD))) * dx*dy
    end

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx, Ny)
    for iy=1:Ny, ix=1:Nx
        if ix == 1 || iy == 1   # highest frequency
            F[ix,iy] = sqrt(Lx*Ly * PSD[ix,iy]) * randn()
        elseif fx[ix] == 0 || fy[iy] == 0   # zero frequency
            F[ix,iy] = sqrt(Lx*Ly * PSD[ix,iy]) * randn()
        elseif fx[ix] > 0 || fy[iy] > 0   # positive frequencies
            F[ix,iy] = sqrt(Lx*Ly * PSD[ix,iy]) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for iy=2:Ny, ix=2:Nx
        if fx[ix] < 0 || fy[iy] < 0   # negative frequencies
            F[ix,iy] = conj(F[Nx - ix + 2, Ny - iy + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / (dx*dy)

    isodd(Nxin) ? R = R[1:end-1,:] : nothing
    isodd(Nyin) ? R = R[:,1:end-1] : nothing

    check_imag(R)

    return real.(R)
end


function rough(xin, yin, zin; ac=:gauss, sigma=nothing, xi=nothing, seed=nothing)
    Nxin, Nyin, Nzin = length(xin), length(yin), length(zin)
    isodd(Nxin) ? x = evenize(xin) : x = xin
    isodd(Nyin) ? y = evenize(yin) : y = yin
    isodd(Nzin) ? z = evenize(zin) : z = zin

    Nx, Ny, Nz = length(x), length(y), length(z)
    dx, dy, dz = x[2]-x[1], y[2]-y[1], z[2]-z[1]
    Lx, Ly, Lz = x[end]-x[1], y[end]-y[1], z[end]-z[1]

    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))
    fy = FFTW.ifftshift(FFTW.fftfreq(Ny, 1/dy))
    fz = FFTW.ifftshift(FFTW.fftfreq(Nz, 1/dz))

    if ac in (:exp, :gauss)
        if isnothing(sigma) || isnothing(xi)
            error("For '$ac' autocorrelation you have to specify 'sigma' and 'xi'")
        end
        if ac == :exp
            psd = psd_exp
        elseif ac == :gauss
            psd = psd_gauss
        end
        if typeof(xi) <: Real
            xix = xiy = xiz = xi
        else
            xix, xiy, xiz = xi
        end
        PSD = [psd(fxn, fyn, fzn, sigma, xix, xiy, xiz) for fxn=fx, fyn=fy, fzn=fz]
    elseif typeof(ac) <: Function
        PSD = [ac(xn, yn, zn) for xn=x, yn=y, zn=z]
        PSD = FFTW.ifftshift(FFTW.fft(FFTW.fftshift(PSD))) * (dx*dy*dz)
    end

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx, Ny, Nz)
    for iz=1:Nz, iy=1:Ny, ix=1:Nx
        if ix == 1 || iy == 1 || iz == 1   # highest frequency
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD[ix,iy,iz]) * randn()
        elseif fx[ix] == 0 || fy[iy] == 0 || fz[iz] == 0   # zero frequency
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD[ix,iy,iz]) * randn()
        elseif fx[ix] > 0 || fy[iy] > 0 || fz[iz] > 0   # positive frequencies
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD[ix,iy,iz]) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for iz=2:Nz, iy=2:Ny, ix=2:Nx
        if fx[ix] < 0 || fy[iy] < 0 || fz[iz] < 0   # negative frequencies
            F[ix,iy,iz] = conj(F[Nx - ix + 2, Ny - iy + 2, Nz - iz + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / (dx*dy*dz)

    isodd(Nxin) ? R = R[1:end-1,:,:] : nothing
    isodd(Nyin) ? R = R[:,1:end-1,:] : nothing
    isodd(Nzin) ? R = R[:,:,1:end-1] : nothing

    check_imag(R)

    return real.(R)
end


# ------------------------------------------------------------------------------------------
"""
    rough(fname; fov, sigma)

Reads external image file and turns it to the map of heights (brighter is higher)

# Arguments
- `fname::String`: input file name

# Keywords
- `fov::Tuple`: field of view, i.e. tuple with the image width and height in SI units
- `height::Real`: the hight (from absolute minimum to absolute maximum) in SI units
"""
function rough(
    fname::String;
    fov, height, lmargin=0, rmargin=0, bmargin=0, tmargin=0, mavg=1, zero_mean=false,
)
    img = Images.load(fname)
    R = @. color2float(img)

    # Properly orient:
    R = transpose(R)
    R = reverse(R; dims=2)

    # Cut margins:
    Nx, Ny = size(R)
    @assert lmargin + rmargin < Nx
    @assert bmargin + tmargin < Ny
    R = R[begin+lmargin:end-rmargin,begin+bmargin:end-tmargin]

    # Smooth:
    if mavg > 1
        R = moving_average(R, mavg)
    end

    # Normalize to a given height:
    R .= R .- minimum(R)
    R .= R ./ maximum(R)
    @. R = R * height

    if zero_mean
        R .= R .- sum(R)/length(R)
    end

    # Create grid:
    Nx, Ny = size(R)
    Lx, Ly = fov
    x = range(-Lx/2, Lx/2, Nx)
    y = range(-Ly/2, Ly/2, Ny)

    return x, y, R
end


function rough(
    x, y, fname::String;
    fov, height, lmargin=0, rmargin=0, bmargin=0, tmargin=0, mavg=1, zero_mean=false,
)
    xex, yex, Rex = rough(
        fname; fov, height, lmargin, rmargin, bmargin, tmargin, mavg, zero_mean,
    )
    itp = linear_interpolation((xex,yex), Rex)
    return [itp(xi, yi) for xi=x, yi=y]
end


# ******************************************************************************************
function geometry_mask(x, z, R::Function; zdirection=1, zshift=0)
    R = [R(xi) for xi=x]
    return geometry_mask(x, z, R; zdirection, zshift)
end

"""
Generates the 2D geomtry mask for a given surface roughness R(x)
"""
function geometry_mask(x, z, R; zdirection=1, zshift=0)
    if zdirection == 1
        gfunc = (x,z,R) -> z <= R + zshift
    elseif zdirection == -1
        gfunc = (x,z,R) -> z >= -R + zshift
    end
    return [gfunc(x[i], zi, R[i]) for i in eachindex(x), zi=z]
end


function geometry_mask(x, y, z, R::Function; zdirection=1, zshift=0)
    R = [R(xi, yi) for xi=x, yi=y]
    return geometry_mask(x, y, z, R; zdirection, zshift)
end

"""
Generates the 3D geomtry mask for a given surface roughness R(x,y)
"""
function geometry_mask(x, y, z, R; zdirection=1, zshift=0)
    @kernel function geometry_mask_kernel!(gmask, x, y, z, R, gfunc)
        ix, iy, iz = @index(Global, NTuple)
        @inbounds begin
            gmask[ix,iy,iz] = gfunc(x[ix], y[iy], z[iz], R[ix,iy])
        end
    end

    if zdirection == 1
        gfunc = (x,y,z,R) -> z <= R + zshift
    elseif zdirection == -1
        gfunc = (x,y,z,R) -> z >= -R + zshift
    end
    gmask = similar(R, length(x), length(y), length(z))
    backend = get_backend(gmask)
    ndrange = size(gmask)
    geometry_mask_kernel!(backend)(gmask, x, y, z, R, gfunc; ndrange)
    return Array{Int}(collect(gmask))
end


# ******************************************************************************************
# Utilities
# ******************************************************************************************
color2float(color) = Float64(color)
color2float(color::Images.ColorTypes.RGB) = Float64(color.r)
color2float(color::Images.ColorTypes.RGBA) = Float64(color.r)


"""
Extends x by one step to make length(x) be even
"""
function evenize(x)
    Nx = length(x)
    dx = x[2] - x[1]
    return range(x[begin], x[end]+dx, Nx+1)
end


"""
Checks that the imaginary part of R is less that its real part by at least lvl times
"""
function check_imag(R; lvl=1e-3)
    maxRre = maximum(abs, extrema(real, R))
    maxRim = maximum(abs, extrema(imag, R))
    if maxRim > lvl * maxRre
        @warn "Imaginary part of R is too big \n" *
              "maximum(real,R)=$maxRre, maximum(imag,R)=$maxRim"
    end
    return nothing
end


"""
Smooths array A using moving average with m neghbors

https://julialang.org/blog/2016/02/iteration/#a_multidimensional_boxcar_filter
"""
function moving_average(A::AbstractArray, m::Int)
    if eltype(A) == Int
        out = zeros(size(A))
    else
        out = similar(A)
    end
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = div(m,2) * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end


end
