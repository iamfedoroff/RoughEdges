import RoughEdges


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 50e-9   # (m) standard deviation (average height)
xi = 50e-9   # (m) correlation length

xmin, xmax, Nx = -1e-6, 1e-6, 1001


# ******************************************************************************************
x = range(xmin, xmax, Nx)

R = RoughEdges.rough(x; sigma, xi, seed)
Ravg = sum(R) / length(R)
Rrms = RoughEdges.rms(R)
@show Ravg/1e-9
@show Rrms/1e-9


# ******************************************************************************************
xu = 1e-6
Ru = 1e-9

RoughEdges.rough_plot(x, R; xu, Ru, new_window=false)
