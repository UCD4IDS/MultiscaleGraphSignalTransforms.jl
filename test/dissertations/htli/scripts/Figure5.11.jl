cd(@__DIR__); include("setups/simpletree.jl");
gr(dpi=200)

for l = 97:100
    plot(size = (200, 500), framestyle = :none)
    gplot!(W, X, width = 1)
    scatter_gplot!(X; marker = 𝚽[:, l])
    plt = plot!(cbar = false, clims = (-0.3, 0.3))
    savefig(plt, "../figs/SimpleTree_eigenvector$(l-1)_red.png")
end
