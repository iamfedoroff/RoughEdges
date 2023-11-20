using RoughEdges
import Statistics: mean, std


# ******************************************************************************************
seed = 1   # random numbers seed

sigma = 1   # (m) standard deviation (average height)
xi = 0.5   # (m) correlation length

# ac(x) = sigma^2 * exp(-abs(x/xi))
ac(x) = sigma^2 * exp(-(x/xi)^2)

xmin, xmax, Nx = -10, 10, 1001


# ******************************************************************************************
x = range(xmin, xmax, Nx)

R = rough(x; sigma, xi, seed)
# R = rough(x; ac=:exp, sigma, xi, seed)
# R = rough(x; ac, seed)

@show mean(R)
@show std(R)


# ******************************************************************************************
rough_plot(x, R)
