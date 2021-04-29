cd(@__DIR__); include("setups/rgc100.jl");
pyplot(dpi = 200)

##
ϵ = 0.3
for J = 1:2
    Uf = unitary_folding_operator(W, GP; ϵ = ϵ, J = J)
    gplot(W, X; width = 1)
    scatter_gplot!(X; marker = diag(Uf), plotOrder = :l2s, ms = 4)
    plt = plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)",
          clims = (0.6, 1), frame = :box, size = (520, 500))
    savefig(plt, "../figs/RGC100_unitary_folding_diag_J$(J).png")
end
