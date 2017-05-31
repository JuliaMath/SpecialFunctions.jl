# Example of Bernoulli functions
#   Replicate (approximately) the plot from http://mathworld.wolfram.com/BernoulliPolynomial.html

using Plots; plotlyjs()
default(show = true) # shouldn't need display commands

using SpecialFunctions
const SF = SpecialFunctions

x = -0.4:0.01:1.0
plt = plot( x, SF.bernoulli.(1,x), 
           label="B_1(x)", legendfont=font(20, "Courier"), legend=:bottomright,
           xaxis = ("x", font(20, "Courier")),
           yaxis = ("", font(20, "Courier"),  (-0.15, 0.15), -0.15:0.05:1.5),
           size = (1200, 800)
           )
 
for i=2:5
    plot!( x, SF.bernoulli.(i,x), label="B_$i(x)")
end

savefig("bernoulli_ex.pdf")
