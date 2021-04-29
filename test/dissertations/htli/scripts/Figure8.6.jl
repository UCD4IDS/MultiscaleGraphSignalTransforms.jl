cd(@__DIR__); include("setups/rgc100.jl");

##
ϵ = 0.3
pair_inds, 𝛟1 = find_pairinds(W; ϵ = ϵ)


## (a) fiedler vector on RGC100
gr(dpi = 200)
plt = gplot(W, X; width = 1)
    scatter_gplot!(X; marker = 𝛟1, plotOrder = :l2s)
    plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)", frame = :box,
    size = (520, 500), right_margin = 5mm)
savefig(plt, "../figs/RGC100_fiedler.png")


## (b) fiedler embedding
pyplot(dpi = 200)
plt = plot(size = (520, 150))
    plot!(-0.055:0.11:0.055, zeros(2), c = :black, lw = 2, yticks = false,
    xlab = latexstring("\\phi^{rw}_1(v)"))
    scatter!(𝛟1, 0:0, ylim = [-0.2, 0.2], grid = false,
    ms = 4, mswidth = 0, legend = false, frame = :box, c = :black)
savefig(plt, "../figs/RGC100_fiedler_embedding.png")
