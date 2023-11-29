using RoughEdges
import Statistics: mean, std


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 1   # (m) standard deviation (average height)
xix = 0.5   # (m) correlation length along x
xiy = 0.5   # (m) correlation length along y

# ac(x, y) = sigma^2 * exp(-abs(x/xix)) * exp(-abs(y/xiy))
ac(x, y) = sigma^2 * exp(-(x/xix)^2) * exp(-(y/xiy)^2)

xmin, xmax, Nx = -10, 10, 501
ymin, ymax, Ny = -10, 10, 501


# ******************************************************************************************
x = range(xmin, xmax, Nx)
y = range(ymin, ymax, Ny)

R = rough(x, y; sigma, xix, xiy, seed)
# R = rough(x, y; ac=:exp, sigma, xix, xiy, seed)
# R = rough(x, y; ac, seed)

@show mean(R)
@show std(R)


# ******************************************************************************************
z = range(-5, 5, 251)

gmask = geometry_mask(x, y, z, R)


# ******************************************************************************************
# using MaxwellPlots

# rough_plot(x, y, R; colorrange=(-4*sigma,4*sigma))

# plot_geometry(x, y, z, gmask)
