cd(@__DIR__); runapprox = false; include("setups/sunflower.jl")
gr(dpi = 200)

## (a) sunflower graph
gplot(W, X; width = 1)
scatter_gplot!(X; c = :red, ms = LinRange(1, 9, N))
plt = plot!(frame = :none)
savefig(plt, "../figs/SunFlower.png")

## (b) voronoi tessellation
xx, yy = getplotxy(voronoiedges(tess))
plt = plot(xx, yy, xlim = [1, 2], ylim = [1, 2], linestyle = :auto, linewidth = 1,
    linecolor = :blue, grid = false, label = "", ratio = 1, frame = :box, ticks = false)
savefig(plt, "../figs/Sunflower_Barbara_Voronoi_cells.png")
