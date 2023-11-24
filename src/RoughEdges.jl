module RoughEdges

import FFTW
import GLMakie as mak
import Images
import KernelAbstractions: @index, @kernel, get_backend
import Random

export rough, rough_plot, geometry_mask


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


function rough(xin, yin; ac=:gauss, sigma=nothing, xix=nothing, xiy=nothing, seed=nothing)
    Nxin, Nyin = length(xin), length(yin)
    isodd(Nxin) ? x = evenize(xin) : x = xin
    isodd(Nyin) ? y = evenize(yin) : y = yin

    Nx, Ny = length(x), length(y)
    dx, dy = x[2]-x[1], y[2]-y[1]
    Lx, Ly = x[end]-x[1], y[end]-y[1]

    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))
    fy = FFTW.ifftshift(FFTW.fftfreq(Ny, 1/dy))

    if ac in (:exp, :gauss)
        if isnothing(sigma) || isnothing(xix) || isnothing(xiy)
            error("For '$ac' autocorrelation you have to specify 'sigma', 'xix' and 'xiy'")
        end
        if ac == :exp
            psd = psd_exp
        elseif ac == :gauss
            psd = psd_gauss
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


function rough(
    xin, yin, zin;
    ac=:gauss, sigma=nothing, xix=nothing, xiy=nothing, xiz=nothing, seed=nothing,
)
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
        if isnothing(sigma) || isnothing(xix) || isnothing(xiy) || isnothing(xiz)
            error("For '$ac' autocorrelation you have to specify 'sigma', 'xix', 'xiy' and 'xiz'")
        end
        if ac == :exp
            psd = psd_exp
        elseif ac == :gauss
            psd = psd_gauss
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
# Visualize rough edges, surfaces, and volumes
# ******************************************************************************************
function rough_plot(x, R; xu=1, Ru=1, new_window=false)
    fig = mak.Figure()
    ax = mak.Axis(fig[1,1]; xlabel="x", ylabel="R")
    mak.lines!(ax, x/xu, R/Ru)
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function rough_plot(
    x, y, R; xu=1, yu=1, Ru=1, colormap=:seismic, colorrange=nothing, new_window=false,
)
    Nx, Ny = size(R)

    if isnothing(colorrange)
        Rmax = maximum(abs, R/Ru)
        colorrange = (-Rmax, Rmax)
    end

    x1, x2 = x[1]/xu, x[end]/xu
    y1, y2 = y[1]/yu, y[end]/yu

    Rx = R[:,div(Ny,2)] / Ru
    Ry = R[div(Nx,2),:] / Ru

    fig = mak.Figure(size=(950,992))
    ax11 = mak.Axis(fig[1,1])
    ax21 = mak.Axis(fig[2,1]; xlabel="x", ylabel="y")
    ax22 = mak.Axis(fig[2,2])

    mak.xlims!(ax11, (x1,x2))
    mak.ylims!(ax11, colorrange...)
    mak.xlims!(ax22, colorrange...)
    mak.ylims!(ax22, (y1,y2))

    mak.linkxaxes!(ax11, ax21)
    mak.linkyaxes!(ax22, ax21)
    mak.hidexdecorations!(ax11, ticks=false)
    mak.hideydecorations!(ax22, ticks=false)

    mak.rowsize!(fig.layout, 1, mak.Relative(1/4))
    mak.colsize!(fig.layout, 1, mak.Relative(3/4))

    mak.lines!(ax11, x/xu, real.(Rx))

    hm = mak.heatmap!(ax21, x/xu, y/yu, real.(R)/Ru; colormap, colorrange)
    mak.lines!(ax21, [x1,x2], [0,0]; color=:black, linewidth=1)
    mak.lines!(ax21, [0,0], [y1,y2]; color=:black, linewidth=1)
    mak.Colorbar(fig[3,1], hm; vertical=false, flipaxis=false)

    mak.lines!(ax22, real.(Ry), y/yu)

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    return nothing
end


function rough_plot(
    x, y, z, R;
    xu=1, yu=1, zu=1, Ru=1, colormap=:seismic, colorrange=nothing, new_window=false,
    aspect=:data, colorbar=true, algorithm=:volumeslices, absorption=1, isovalue=0,
)
    Nx, Ny, Nz = size(R)

    if isnothing(colorrange)
        Rmax = maximum(abs, R/Ru)
        colorrange = (-Rmax, Rmax)
    end

    ix, iy, iz = div(Nx,2), div(Ny,2), div(Nz,2)

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis3(fig[1,1]; xlabel="x", ylabel="y", zlabel="z", aspect, perspectiveness=0)

    if algorithm == :volumeslices
        vol = mak.volumeslices!(
            ax, x/xu, y/yu, z/zu, R/Ru; bbox_visible=false, colormap, colorrange,
        )
        vol.update_yz[](ix)
        vol.update_xz[](iy)
        vol.update_xy[](iz)

        if colorbar
            mak.Colorbar(fig[1,2], vol.heatmap_xy[], height=mak.Relative(0.5))
        end

        sg = mak.SliderGrid(
            fig[2,1],
            (label="x", range=1:length(x), startvalue=ix),
            (label="y", range=1:length(y), startvalue=iy),
            (label="z", range=1:length(z), startvalue=iz),
        )
        mak.on(sg.sliders[1].value) do i
            vol.update_yz[](i)
        end
        mak.on(sg.sliders[2].value) do i
            vol.update_xz[](i)
        end
        mak.on(sg.sliders[3].value) do i
            vol.update_xy[](i)
        end

        # hmaps = [vol.heatmap_yz[], vol.heatmap_xz[], vol.heatmap_xy[]]
        # toggles = [mak.Toggle(sg.layout[i,4], active=true) for i in 1:length(hmaps)]
        # map(zip(hmaps, toggles)) do (h, t)
        #     mak.connect!(h.visible, t.active)
        # end
    else
        Rmax = maximum(abs, R/Ru)

        img = mak.volume!(
            ax, x/xu, y/yu, z/zu, R/Ru;
            colormap, colorrange, algorithm, absorption, isovalue, isorange=0.05*Rmax,
        )

        if colorbar
            mak.Colorbar(fig[1,2], img, height=mak.Relative(0.5))
        end
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    return nothing
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
