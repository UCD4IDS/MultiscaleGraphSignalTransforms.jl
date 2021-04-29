cd(@__DIR__); include("setups/rgc100.jl")
plotlyjs(dpi = 200)

gplot(W, X3; grid = true, shape = :circle, mcolor = :yellow)
    plot!(xlab = "x", ylab = "y", zlab = "z", size = (500, 500))
