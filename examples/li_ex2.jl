# do plots of the Li function,
#  this one takes a while

using Plots; plotlyjs()
using PlotUtils  # already from Plots, but as a reminder
using PlotThemes # already from Plots, but as a reminder
default(show = true) # shouldn't need display commands

using SpecialFunctions
const SF = SpecialFunctions

# one of the figures from https://en.wikipedia.org/wiki/Polylogarithm
#   presuming this is a phase plot (there is not scale)
# step = 0.02   $ took about 200 seconds
# step = 0.01   # took about 11 minutes
step = 0.0025 
xs = -2.0 : step : 2.0
ys = -2.0 : step : 2.0
Z = [Complex(x, y) for x in xs, y in ys]'
tic()
L = SF.Li.(-2,Z)
theTime = toc()
a = angle(L)/pi
i = abs(imag(L)) .< 1.0e-8
a[i] = (-sign(real(L[i]))+1.0)/2.0
j = abs(L) .< 1.0e-8
a[j] = NaN
heatmap(xs, ys, a,
        xaxis = ("real", font(20, "Courier")),
        yaxis = ("imag", font(20, "Courier")),
        size = (800, 800),
        color = :Spectral)

sleep(1)
savefig("li_ex2.pdf")
