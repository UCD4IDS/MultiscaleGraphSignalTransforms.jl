cd(@__DIR__); include("setups/path16.jl");
gr(dpi = 250)

plot(size = (800, 250), framestyle = :none, xlim = [-10.5, N+13], ylim = [-5, 1])
path_pc_plot!(W, X; c = :yellow)

k = 1
path_pc_plot!(W, X .- [12 4]; c = :yellow,
        ind = GP_primal.inds[GP_primal.rs[k, 2]:GP_primal.rs[k+1, 2] - 1, 2])

k = 2
path_pc_plot!(W, X .- [-12 4]; c = :yellow,
        ind = GP_primal.inds[GP_primal.rs[k, 2]:GP_primal.rs[k+1, 2] - 1, 2])

plot!([6, -1], [-0.8, -3.3], line = :solid, c = :black, lw = 0.65)
plot!([11, 18], [-0.8, -3.3], line = :solid, c = :black, lw = 0.65)
annotate!(mean(X[:, 1]), 0.5,
    text(latexstring("V= \\{ \\delta_{1}, ..., \\delta_{16} \\}"), :bottom, 11))
annotate!(mean(X[:, 1]) - 12, -4 - 0.7,
    text(latexstring("V_{0}= \\{ \\delta_{1}, \\delta_{3}, \\delta_{5}, \\delta_{7}, \\delta_{10}, \\delta_{12}, \\delta_{14}, \\delta_{16} \\}"), :top, 11))
annotate!(mean(X[:, 1]) + 12, -4 - 0.7,
    text(latexstring("V_{1}= \\{ \\delta_{2}, \\delta_{4}, \\delta_{6}, \\delta_{8}, \\delta_{9}, \\delta_{11}, \\delta_{13}, \\delta_{15} \\}"), :top, 11))

plt = current()
savefig(plt, "../figs/Path16_1lev_PC.png")
display(plt)
