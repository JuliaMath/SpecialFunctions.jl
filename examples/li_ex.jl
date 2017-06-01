# do plots of the Li function,



using Plots; plotlyjs()
default(show = true) # shouldn't need display commands
theme(:solarized_light)

using SpecialFunctions
const SF = SpecialFunctions

# replicate (approximately) the plot from http://mathworld.wolfram.com/Polylogarithm.html
x = -2.0 : 0.01 : 1.0
plt = plot( x, real(SF.Li.(-2,x)), 
           label="Li_{-2}(x)", legendfont=font(20, "Courier"), legend=:bottomright,
           xaxis = ("x", font(20, "Courier")),
           yaxis = ("", font(20, "Courier"),  (-1, 1), -1.0:0.25:1.0),
           size = (1200, 800)
           )
 
for i=-1:2
    plot!( x, real(SF.Li.(i,x)), label="Li_{$i}(x)")
end

sleep(1)
savefig("li_ex.pdf")


