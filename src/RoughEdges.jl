module RoughEdges

import FFTW
import GLMakie as mak
import KernelAbstractions: @index, @kernel, get_backend
import Random

export rough, rough_plot, rms


# ******************************************************************************************
# Power spectral density (PSD) functions
#
# C.A. Mack, "Analytic form for the power spectral density in one, two, and three
# dimensions", Journal of Micro/Nanolitography, MEMS, and MOEMS, 10, 040501 (2011)
# ******************************************************************************************
psd1D(f, sigma, xi) = sqrt(pi) * sigma^2 * xi * exp(-(pi * f * xi)^2)

psd2D(f, sigma, xi) = pi * sigma^2 * xi^2 * exp(-(pi * f * xi)^2)

psd3D(f, sigma, xi) = pi^(3/2) * sigma^2 * xi^3 * exp(-(pi * f * xi)^2)


# ******************************************************************************************
# Rough edge, surface, and volume generators
#
# C.A. Mack, "Generating random rough edges, surfaces, and volumes", Applied Optics, 52,
# 1472 (2013); https://doi.org/10.1364/AO.52.001472
# ******************************************************************************************
function rough(xin; sigma, xi, seed=nothing, psd=psd1D)
    Nxin = length(xin)
    isodd(Nxin) ? x = evenize(xin) : x = xin

    Nx = length(x)
    dx = x[2]-x[1]
    Lx = x[end]-x[1]

    # Spectral grid:
    # f[1]   # highest frequency
    # f[2:N/2]   # negative frequencies
    # f[N/2+1]   # zero frequency
    # f[N/2+2:N]   # positive frequencies
    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx)
    for ix=1:Nx
        f = fx[ix]
        PSD = psd(f, sigma, xi)   # power spectral density
        if ix == 1   # highest frequency
            F[ix] = sqrt(Lx * PSD) * randn()
        elseif fx[ix] == 0   # zero frequency
            F[ix] = sqrt(Lx * PSD) * randn()
        elseif fx[ix] > 0   # positive frequencies
            F[ix] = sqrt(Lx * PSD) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for ix=2:Nx
        if fx[ix] < 0   # negative frequencies
            F[ix] = conj(F[Nx - ix + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / dx
    maxRre = maximum(abs, extrema(real, R))
    maxRim = maximum(abs, extrema(imag, R))
    if maxRim > 1e-6*maxRre
        @warn "max(imag(R)) > 1e-6*max(real(R))"
        @show maxRre, maxRim
    end

    isodd(Nxin) ? R = R[1:end-1] : nothing

    return real.(R)
end


function rough(xin, yin; sigma, xi, seed=nothing, psd=psd2D)
    Nxin, Nyin = length(xin), length(yin)
    isodd(Nxin) ? x = evenize(xin) : x = xin
    isodd(Nyin) ? y = evenize(yin) : y = yin

    Nx, Ny = length(x), length(y)
    dx, dy = x[2]-x[1], y[2]-y[1]
    Lx, Ly = x[end]-x[1], y[end]-y[1]

    # Spectral grid:
    # f[1]   # highest frequency
    # f[2:N/2]   # negative frequencies
    # f[N/2+1]   # zero frequency
    # f[N/2+2:N]   # positive frequencies
    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))
    fy = FFTW.ifftshift(FFTW.fftfreq(Ny, 1/dy))

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx, Ny)
    for iy=1:Ny, ix=1:Nx
        f = sqrt(fx[ix]^2 + fy[iy]^2)
        PSD = psd(f, sigma, xi)   # power spectral density
        if ix == 1 || iy == 1   # highest frequency
            F[ix,iy] = sqrt(Lx*Ly * PSD) * randn()
        elseif fx[ix] == 0 || fy[iy] == 0   # zero frequency
            F[ix,iy] = sqrt(Lx*Ly * PSD) * randn()
        elseif fx[ix] > 0 || fy[iy] > 0   # positive frequencies
            F[ix,iy] = sqrt(Lx*Ly * PSD) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for iy=2:Ny, ix=2:Nx
        if fx[ix] < 0 || fy[iy] < 0   # negative frequencies
            F[ix,iy] = conj(F[Nx - ix + 2, Ny - iy + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / (dx*dy)
    maxRre = maximum(abs, extrema(real, R))
    maxRim = maximum(abs, extrema(imag, R))
    if maxRim > 1e-6*maxRre
        @warn "max(imag(R)) > 1e-6*max(real(R))"
        @show maxRre, maxRim
    end

    isodd(Nxin) ? R = R[1:end-1,:] : nothing
    isodd(Nyin) ? R = R[:,1:end-1] : nothing

    return real.(R)
end


function rough(xin, yin, zin; sigma, xi, seed=nothing, psd=psd3D)
    Nxin, Nyin, Nzin = length(xin), length(yin), length(zin)
    isodd(Nxin) ? x = evenize(xin) : x = xin
    isodd(Nyin) ? y = evenize(yin) : y = yin
    isodd(Nzin) ? z = evenize(zin) : z = zin

    Nx, Ny, Nz = length(x), length(y), length(z)
    dx, dy, dz = x[2]-x[1], y[2]-y[1], z[2]-z[1]
    Lx, Ly, Lz = x[end]-x[1], y[end]-y[1], z[end]-z[1]

    # Spectral grid:
    # f[1]   # highest frequency
    # f[2:N/2]   # negative frequencies
    # f[N/2+1]   # zero frequency
    # f[N/2+2:N]   # positive frequencies
    fx = FFTW.ifftshift(FFTW.fftfreq(Nx, 1/dx))
    fy = FFTW.ifftshift(FFTW.fftfreq(Ny, 1/dy))
    fz = FFTW.ifftshift(FFTW.fftfreq(Nz, 1/dz))

    isnothing(seed) ? nothing : Random.seed!(seed)

    F = zeros(ComplexF64, Nx, Ny, Nz)
    for iz=1:Nz, iy=1:Ny, ix=1:Nx
        f = sqrt(fx[ix]^2 + fy[iy]^2 + fz[iz]^2)
        PSD = psd(f, sigma, xi)   # power spectral density
        if ix == 1 || iy == 1 || iz == 1   # highest frequency
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD) * randn()
        elseif fx[ix] == 0 || fy[iy] == 0 || fz[iz] == 0   # zero frequency
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD) * randn()
        elseif fx[ix] > 0 || fy[iy] > 0 || fz[iz] > 0   # positive frequencies
            F[ix,iy,iz] = sqrt(Lx*Ly*Lz * PSD) * (randn() + 1im*randn()) / sqrt(2)
        end
    end
    for iz=2:Nz, iy=2:Ny, ix=2:Nx
        if fx[ix] < 0 || fy[iy] < 0 || fz[iz] < 0   # negative frequencies
            F[ix,iy,iz] = conj(F[Nx - ix + 2, Ny - iy + 2, Nz - iz + 2])
        end
    end

    R = FFTW.ifftshift(FFTW.ifft(FFTW.fftshift(F))) / (dx*dy*dz)
    maxRre = maximum(abs, extrema(real, R))
    maxRim = maximum(abs, extrema(imag, R))
    if maxRim > 1e-6*maxRre
        @warn "max(imag(R)) > 1e-6*max(real(R))"
        @show maxRre, maxRim
    end

    isodd(Nxin) ? R = R[1:end-1,:,:] : nothing
    isodd(Nyin) ? R = R[:,1:end-1,:] : nothing
    isodd(Nzin) ? R = R[:,:,1:end-1] : nothing

    return real.(R)
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

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax11 = mak.Axis(fig[1,1])
    ax21 = mak.Axis(fig[2,1]; xlabel="x", ylabel="y")
    ax22 = mak.Axis(fig[2,2])

    mak.ylims!(ax11, colorrange...)
    mak.ylims!(ax22, colorrange...)

    mak.linkxaxes!(ax11, ax21)
    mak.linkyaxes!(ax22, ax21)
    mak.hidexdecorations!(ax11, ticks=false)
    mak.hideydecorations!(ax22, ticks=false)

    mak.rowsize!(fig.layout, 1, mak.Relative(1/4))
    mak.colsize!(fig.layout, 1, mak.Relative(3/4))

    mak.lines!(ax11, x/xu, real.(Rx))

    hm = mak.heatmap!(ax21, x/xu, y/yu, real.(R)/Ru; colormap, colorrange)
    mak.lines!(ax21, [x1,x2], [0,0]; color=:black, linewidth=0.5)
    mak.lines!(ax21, [0,0], [y1,y2]; color=:black, linewidth=0.5)
    mak.Colorbar(fig[3,1], hm; vertical=false, flipaxis=false)

    mak.lines!(ax22, real.(Ry), y/yu)
    mak.ylims!(ax22, (x1,x2))

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

    fig = mak.Figure(resolution=(950,992), fontsize=14)
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
"""
Extends x by one step to make length(x) be even
"""
function evenize(x)
    Nx = length(x)
    dx = x[2] - x[1]
    return range(x[begin], x[end]+dx, Nx+1)
end


"""
Calculates root-mean-square deviation
"""
function rms(x)
    xavg = sum(x) / length(x)
    tmp = zero(eltype(x))
    for i in eachindex(x)
        tmp += (x[i] - xavg)^2
    end
    return sqrt(tmp / length(x))
end


"""
Prepares geomtry mask for a given grid (x,y,z) according to the heights distribution R and
logical function gfunc.

    Example:
    gfunc(x,y,z,R) = (z - 1e-6) >= R   # rough  surface located at z=1um
"""
@kernel function geometry_mask_kernel!(gmask, x, y, z, R, gfunc)
    ix, iy, iz = @index(Global, NTuple)
    @inbounds begin
        gmask[ix,iy,iz] = gfunc(x[ix], y[iy], z[iz], R[ix,iy])
    end
end
function geometry_mask(x, y, z, R; gfunc)
    gmask = similar(R, length(x), length(y), length(z))
    backend = get_backend(gmask)
    ndrange = size(gmask)
    geometry_mask_kernel!(backend)(gmask, x, y, z, R, gfunc; ndrange)
    return Array{Int}(collect(gmask))
end


end