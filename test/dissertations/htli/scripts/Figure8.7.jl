cd(@__DIR__); include("setups/rgc100.jl");
pyplot(dpi = 200)

##
ϵ = 0.3
pair_inds, 𝛟1 = find_pairinds(W; ϵ = ϵ)
v_bar = (𝛟1[pair_inds[1, end-1:end]] - 𝛟1[pair_inds[2, end-1:end]]) ./ 2
band = ϵ * norm(𝛟1, Inf)
pos_ar = findall(0 .< 𝛟1 .< band)
neg_ar = findall(-band .< 𝛟1 .<= 0)

# soft partition 1D
plt = plot(layout = Plots.grid(2, 1), size = (520, 500))
    plot!(-0.055:0.11:0.055, zeros(2), c = :black, lw = 2, yticks = false,
    subplot = 1)
    plot!(zeros(2), [-0.1, 0.1], lw = 3, c = :grey, linestyle = :solid,
    xlab = latexstring("\\phi^{rw}_1(v)"), subplot = 1)
    scatter!(𝛟1[𝛟1 .> 0], 0:0, ylim = [-0.2, 0.2], grid = false,
    ms = 4, mswidth = 0,  legend = false, frame = :box, c = :red, subplot = 1)
    scatter!(𝛟1[𝛟1 .< 0], 0:0, m = (1, 4), mswidth = 0, c = :blue, subplot = 1)
    plot!(ones(2) * ϵ * norm(𝛟1, Inf), [-0.08, 0.08], lw = 2, c = :grey,
    linestyle = :dash, subplot = 1)
    plot!(-ones(2) * ϵ * norm(𝛟1, Inf), [-0.08, 0.08], lw = 2, c = :grey,
    linestyle = :dash, subplot = 1)
    annotate!(0.02, -0.12,
    text(latexstring("\\epsilon \\cdot \\Vert \\phi^{rw}_1 \\Vert_{\\infty}"), 11),
    subplot = 1)
    annotate!(-0.02, -0.12,
    text(latexstring("-\\epsilon \\cdot \\Vert \\phi^{rw}_1 \\Vert_{\\infty}"), 11),
    subplot = 1)
    annotate!(0, 0.15, text("action region", 10), subplot = 1)

    plot!(-0.016:0.032:0.016, zeros(2), c = :black, lw = 2, yticks = false,
    subplot = 2)
    plot!(zeros(2), [-0.1, 0.1], lw = 3, c = :grey, linestyle = :solid,
    xlab = latexstring("\\phi^{rw}_1(v)"), subplot = 2)
    scatter!(𝛟1[pos_ar], 0:0, ylim = [-0.2, 0.2], grid = false,
    ms = 4, mswidth = 0,  legend = false, frame = :box, c = :red, subplot = 2)
    scatter!(𝛟1[neg_ar], 0:0, m = (1, 4), mswidth = 0, c = :blue, subplot = 2)
    plot!(ones(2) * ϵ * norm(𝛟1, Inf), [-0.08, 0.08], lw = 2, c = :grey,
    linestyle = :dash, subplot = 2)
    plot!(-ones(2) * ϵ * norm(𝛟1, Inf), [-0.08, 0.08], lw = 2, c = :grey,
    linestyle = :dash, subplot = 2)
    annotate!(0.012, -0.12,
    text(latexstring("\\epsilon \\cdot \\Vert \\phi^{rw}_1 \\Vert_{\\infty}"), 11),
    subplot = 2)
    annotate!(-0.012, -0.12,
    text(latexstring("-\\epsilon \\cdot \\Vert \\phi^{rw}_1 \\Vert_{\\infty}"), 11),
    subplot = 2)
    annotate!(-0.0075, 0.15, text("negative action region", 10), subplot = 2)
    annotate!(0.0075, 0.15, text("positive action region", 10), subplot = 2)
savefig(plt, "../figs/RGC100_fiedler_embedding_soft_clustering.png")

# zoom-in: find relection points about the nodal point
plt = plot(layout = Plots.grid(2, 1), size = (520, 500))
    plot!(-0.005:0.01:0.005, zeros(2), c = :black, lw = 2, yticks = false,
    subplot = 1)
    scatter!(𝛟1[pair_inds[1, end-1:end]], 0:0, ylim = [-0.2, 0.2], grid = false,
    ms = 6, mswidth = 0,  legend = false, frame = :box, c = :red, subplot = 1)
    scatter!(𝛟1[pair_inds[2, end-1:end]], 0:0, ms = 6, mswidth = 0, c = :blue,
    subplot = 1)
    plot!(zeros(2), [-0.1, 0.1], lw = 3, c = :grey, linestyle = :solid,
    xlab = latexstring("\\phi^{rw}_1(v)"), subplot = 1)
    annotate!(𝛟1[pair_inds[1, end]], 0.05, text(latexstring("v^{+}_{1}"), 11),
    subplot = 1)
    annotate!(𝛟1[pair_inds[1, end-1]], 0.05, text(latexstring("v^{+}_{2}"), 11),
    subplot = 1)
    annotate!(𝛟1[pair_inds[2, end]], 0.05, text(latexstring("v^{-}_{1}"), 11),
    subplot = 1)
    annotate!(𝛟1[pair_inds[2, end-1]], 0.05, text(latexstring("v^{-}_{2}"), 11),
    subplot = 1)

    plot!(vec(-0.005:0.01:0.005), zeros(2), c = :black, lw = 2,
    yticks = false, subplot = 2)
    scatter!(vcat(v_bar, -v_bar), 0:0, ylim = [-0.2, 0.2], grid = false,
    ms = 6, mswidth = 0,  legend = false, frame = :box, c = :purple, subplot = 2)
    plot!(zeros(2), [-0.1, 0.1], lw = 3, c = :grey, linestyle = :solid,
    xlab = latexstring("r"), subplot = 2)
    annotate!(v_bar[end], -0.05, text(latexstring("r_{1}"), 11),
    subplot = 2)
    annotate!(v_bar[end-1], -0.05, text(latexstring("r_{2}"), 11),
    subplot = 2)
    annotate!(-v_bar[end], -0.05, text(latexstring("-r_{1}"), 11),
    subplot = 2)
    annotate!(-v_bar[end-1], -0.05, text(latexstring("-r_{2}"), 11),
    subplot = 2)
savefig(plt, "../figs/RGC100_fiedler_embedding_find_reflection_points.png")
