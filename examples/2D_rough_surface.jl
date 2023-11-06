import RoughEdges


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 50e-9   # (m) standard deviation (average height)
xi = 50e-9   # (m) correlation length

xmin, xmax, Nx = -5e-6, 5e-6, 1001
ymin, ymax, Ny = -5e-6, 5e-6, 1001


# ******************************************************************************************
x = range(xmin, xmax, Nx)
y = range(ymin, ymax, Ny)

R = RoughEdges.rough(x, y; sigma, xi, seed)
Ravg = sum(R) / length(R)
Rrms = RoughEdges.rms(R)
@show Ravg/1e-9
@show Rrms/1e-9


# ******************************************************************************************
xu = yu = 1e-6
Ru = 1e-9

colorrange = (-4*sigma/Ru, 4*sigma/Ru)

RoughEdges.rough_plot(x, y, R; xu, yu, Ru, colorrange, new_window=false)



# ******************************************************************************************
# zmin, zmax, Nz = -1.1e-6, 0.3e-6, 141
# z = range(zmin, zmax, Nz)

# gfunc(x, y, z, R) = z >= R   # gmask function

# gmask = RoughEdges.geometry_mask(x, y, z, R; gfunc)
# @show sum(gmask)

# using MaxwellPlots
# xu = yu = zu = 1e-6
# plot_geometry(x, y, z, gmask; xu, yu, zu)
