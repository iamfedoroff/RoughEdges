import RoughEdges


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 50e-9   # (m) standard deviation (average height)
xi = 50e-9   # (m) correlation length

xmin, xmax, Nx = -1e-6, 1e-6, 251
ymin, ymax, Ny = -1e-6, 1e-6, 251
zmin, zmax, Nz = -1e-6, 1e-6, 251


# ******************************************************************************************
x = range(xmin, xmax, Nx)
y = range(ymin, ymax, Ny)
z = range(zmin, zmax, Nz)

R = RoughEdges.rough(x, y, z; sigma, xi, seed)
Ravg = sum(R) / length(R)
Rrms = RoughEdges.rms(R)
@show Ravg/1e-9
@show Rrms/1e-9


# ******************************************************************************************
xu = yu = zu = 1e-6
Ru = 1e-9

colorrange = (-4*sigma/Ru, 4*sigma/Ru)

RoughEdges.rough_plot(x, y, z, R; xu, yu, zu, Ru, colorrange, new_window=false)
