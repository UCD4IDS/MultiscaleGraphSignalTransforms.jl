cd(@__DIR__); include("setups/path128.jl");
gr(dpi = 200)

## (a)
plt = plot(size = (400, 480), layout = Plots.grid(2, 1), ylim = [-0.3, 1.3])
    scatter!(χ(1:64, N), c = :black, grid = false, legend = false, frame = :box,
    lw = 1.5, ms = 3)
    bar!(χ(1:64, N); bar_width = 0.02, c = :black)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))

    scatter!(χ(65:N, N), c = :black, grid = false, legend = false, frame = :box,
    lw = 1.5, subplot = 2, ms = 3)
    bar!(χ(65:N, N); bar_width = 0.02, c = :black, subplot = 2)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]),
    subplot = 2)
savefig(plt, "../figs/Path128_hard_bipart_restric_j1.png")

## (b)
plt = plot(size = (400, 480), layout = Plots.grid(2, 1), ylim = [-0.3, 1.3])
    scatter!(P0 * ones(N), c = :black, grid = false, legend = false, frame = :box,
    lw = 1.5, ms = 3)
    bar!(P0 * ones(N); bar_width = 0.02, c = :black)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))

    scatter!(P1 * ones(N), c = :black, grid = false, legend = false, frame = :box,
    lw = 1.5, ms = 3, subplot = 2)
    bar!(P1 * ones(N); bar_width = 0.02, c = :black, subplot = 2)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]),
    subplot = 2)
savefig(plt, "../figs/Path128_soft_bipart_PSO_j1.png")
