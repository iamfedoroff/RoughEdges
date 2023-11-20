using RoughEdges
import Statistics: mean, std


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 1   # (m) standard deviation (average height)
xix = 0.5   # (m) correlation length along x
xiy = 0.5   # (m) correlation length along y
xiz = 0.5   # (m) correlation length along z

# ac(x, y, z) = sigma^2 * exp(-abs(x/xix)) * exp(-abs(y/xiy)) * exp(-abs(z/xiz))
ac(x, y, z) = sigma^2 * exp(-(x/xix)^2) * exp(-(y/xiy)^2) * exp(-(z/xiz)^2)

xmin, xmax, Nx = -10, 10, 251
ymin, ymax, Ny = -10, 10, 251
zmin, zmax, Nz = -10, 10, 251


# ******************************************************************************************
x = range(xmin, xmax, Nx)
y = range(ymin, ymax, Ny)
z = range(zmin, zmax, Nz)

R = rough(x, y, z; sigma, xix, xiy, xiz, seed)
# R = rough(x, y, z; ac=:exp, sigma, xix, xiy, xiz, seed)
# R = rough(x, y, z; ac, seed)

@show mean(R)
@show std(R)


# ******************************************************************************************
rough_plot(x, y, z, R; colorrange=(-4*sigma,4*sigma))
