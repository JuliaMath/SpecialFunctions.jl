# do plots of the Li function, similar to those of 
#            Vepstas, L. (2008). "An efficient algorithm for accelerating the convergence of oscillatory series, useful for computing the polylogarithm and Hurwitz zeta functions". Numerical Algorithms. 47 (3): 211–252. arXiv:math.CA/0702243Freely accessible. doi:10.1007/s11075-007-9153-8.

using PyPlot
include("li.jl")

# Figure 5
s = 0.5 + 15*im
step = 0.005
x = -3.5: step :3.5
y = -3.5: step :3.5

L = zeros(Complex, (length(x), length(y)) )
for i=1:length(x)
    for j=1:length(y)
        # println("z = $z")
        z = x[i]+im*y[j]
        L[i,j] = Li(s,z)
    end
end

# T = abs(L)
# T[find(abs(T) .> 1000)] = 1000

T = angle(L)

surf(x, y, T, facecolors=get_cmap("jet").o(T/maximum(T)))
xlabel("Re")
ylabel("Im")
