import RoughEdges


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 50e-9   # (m) standard deviation (average height)
xix = 50e-9   # (m) correlation length along x
xiy = 50e-9   # (m) correlation length along y
xiz = 50e-9   # (m) correlation length along z

# ac(x, y, z) = sigma^2 * exp(-abs(x/xix)) * exp(-abs(y/xiy)) * exp(-abs(z/xiz))
ac(x, y, z) = sigma^2 * exp(-(x/xix)^2) * exp(-(y/xiy)^2) * exp(-(z/xiz)^2)

xmin, xmax, Nx = -1e-6, 1e-6, 251
ymin, ymax, Ny = -1e-6, 1e-6, 251
zmin, zmax, Nz = -1e-6, 1e-6, 251


# ******************************************************************************************
x = range(xmin, xmax, Nx)
y = range(ymin, ymax, Ny)
z = range(zmin, zmax, Nz)

R = RoughEdges.rough(x, y, z; sigma, xix, xiy, xiz, seed)
# R = RoughEdges.rough(x, y, z; ac=:exp, sigma, xix, xiy, xiz, seed)
# R = RoughEdges.rough(x, y, z; ac, seed)
Ravg = sum(R) / length(R)
Rrms = RoughEdges.rms(R)
@show Ravg/1e-9
@show Rrms/1e-9


# ******************************************************************************************
xu = yu = zu = 1e-6
Ru = 1e-9

colorrange = (-4*sigma/Ru, 4*sigma/Ru)

RoughEdges.rough_plot(x, y, z, R; xu, yu, zu, Ru, colorrange, new_window=false)
