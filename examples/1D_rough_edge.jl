import RoughEdges


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 50e-9   # (m) standard deviation (average height)
xi = 50e-9   # (m) correlation length

xmin, xmax, Nx = -5e-6, 5e-6, 1001
zmin, zmax, Nz = -1.1e-6, 0.3e-6, 141


# ******************************************************************************************
x = range(xmin, xmax, Nx)
z = range(zmin, zmax, Nz)

R = RoughEdges.rough(x; sigma, xi, seed)
Ravg = sum(R) / length(R)
Rrms = RoughEdges.rms(R)
@show Ravg/1e-9
@show Rrms/1e-9


# ******************************************************************************************
xu = zu = 1e-6
Ru = 1e-9

RoughEdges.rough_plot(x, R; xu, Ru, new_window=false)
