cd(@__DIR__); include("setups/path128.jl");
gr(dpi = 200)

## (a)
plt = plot(size = (500, 400))
    heatmap!(Uf, c = :viridis, ratio = 1, yflip = true, xlim = [1, N], clims = (-1, 1))
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
    yticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
savefig(plt, "../figs/Path128_Uf_j1.png")

## (b)
plt = plot(size = (500, 560), layout = Plots.grid(2, 1))
    scatter!(diag(Uf), ylim = [-1, 1], c = :black, grid = false, legend = false,
    frame = :box, lw = 1.5, ms = 3)
    bar!(diag(Uf); bar_width = 0.02, c = :black)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
    title!("Diagonal entries")

    scatter!(anti_diag(Uf), ylim = [-1, 1], c = :black, grid = false, legend = false,
    frame = :box, lw = 1.5, ms = 3, subplot = 2)
    bar!(anti_diag(Uf); bar_width = 0.02, c = :black, subplot = 2)
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
    title!("Antidiagonal entries", subplot = 2)
savefig(plt, "../figs/Path128_Uf_diag_entries.png")
