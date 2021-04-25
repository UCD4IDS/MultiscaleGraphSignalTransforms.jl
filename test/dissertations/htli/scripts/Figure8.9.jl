cd(@__DIR__); include("setups/rgc100.jl")
gr(dpi = 200)

## construct full HGLET and Lapped-HGLET dictionaries
@time HGLET_dic = HGLET_dictionary(GP, G_Sig; method = :L)
@time LPHGLET_dic = LPHGLET_dictionary(GP, G_Sig; method = :L, ϵ = 0.3)

## generate figures
# pre-selected (j, k, l)s
T = [(2, 1, 1), (2, 1, 12), (3, 1, 3), (3, 1, 5),
    (4, 5, 6), (6, 23, 3), (7, 5, 2), (8, 15, 1)]

for (j, k, l) in T
    𝚽H = HGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'
    𝚽lH = LPHGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'
    p1 = gplot(W, X; width = 1)
        scatter_gplot!(X; marker = 𝚽H[:, l], plotOrder = :s2l)
        plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)", frame = :box, cbar = false)
        plot!(left_margin = 3mm, bottom_margin = 3mm, clims = (-0.05, 0.05))
    p2 = gplot(W, X; width = 1)
        scatter_gplot!(X; marker = 𝚽lH[:, l], plotOrder = :s2l)
        plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)", frame = :box, cbar = false)
        plot!(left_margin = 3mm, bottom_margin = 3mm, clims = (-0.05, 0.05))
    plt = plot(p1, p2, size = (800, 400), layout = Plots.grid(1, 2))
    savefig(plt, "../figs/RGC100_HGLET_LPHGLET_j$(j-1)k$(k-1)l$(l-1).png")
end
