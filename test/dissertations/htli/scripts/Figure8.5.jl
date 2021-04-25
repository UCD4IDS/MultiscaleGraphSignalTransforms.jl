cd(@__DIR__); include("setups/path128.jl");
gr(dpi = 200)

## (a)
Uf = unitary_folding_operator(W, GP; ϵ = 0.38, J = 2)' # ϵa ≈ 0.38⋅∥𝛟₁∥_{∞}
plt = plot(size = (500, 400))
    heatmap!(Uf, c = :viridis, ratio = 1, yflip = true, xlim = [1, N], clims = (-1, 1))
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
    yticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
savefig(plt, "../figs/Path128_Uf_j2.png")

## (b)
Uf = unitary_folding_operator(W, GP; ϵ = 0.38, J = 3)'
plt = plot(size = (500, 400))
    heatmap!(Uf, c = :viridis, ratio = 1, yflip = true, xlim = [1, N], clims = (-1, 1))
    xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
    yticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
savefig(plt, "../figs/Path128_Uf_j3.png")
