cd(@__DIR__); include("setups/rgc100.jl")
gr(dpi = 200)

## construct full HGLET and Lapped-HGLET dictionaries
@time HGLET_dic = HGLET_dictionary(GP, G_Sig; gltype = :L)
@time LPHGLET_dic = LPHGLET_dictionary(GP, G_Sig; gltype = :L, ϵ = 0.3)

## generate figures
# pre-selected (j, k, l)s
T = [(2, 1, 1), (2, 1, 12), (3, 1, 3), (3, 1, 5),
    (4, 5, 6), (6, 23, 3), (7, 5, 2), (8, 15, 1)]

XLs = [
    [-50, 50], [-50, 50], [-110, 40], [-110, 40],
    [-20, 50], [50, 120], [-80, -20], [-40, -10]
]
YLs = [
    [-50, 50], [-50, 50], [-50, 100], [-50, 100],
    [-50, 20], [-70, 0], [-90, -30], [-160, -130]
]

for i in 1:length(T)
    (j, k, l) = T[i]
    𝚽H = HGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'
    𝚽lH = LPHGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'
    p1 = gplot(W, X; width = 1)
        scatter_gplot!(X; marker = 𝚽H[:, l], plotOrder = :s2l, ms = 5)
        plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)", frame = :box, cbar = false)
        plot!(left_margin = 3mm, bottom_margin = 3mm, clims = (-0.05, 0.05))
    p2 = gplot(W, X; width = 1)
        scatter_gplot!(X; marker = 𝚽lH[:, l], plotOrder = :s2l, ms = 5)
        plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)", frame = :box, cbar = false)
        plot!(left_margin = 3mm, bottom_margin = 3mm, clims = (-0.05, 0.05))
    plt = plot(p1, p2, size = (640, 300), layout = Plots.grid(1, 2),
        xlims = XLs[i], ylims = YLs[i])
        # xlims = [-120, 120], ylims = [-150, 150])
    # gplot(W, X; width = 1)
    # ord = :s2l 
    # if i == 8
    #     ord = :l2s
    # end
    # scatter_gplot!(X; marker = 𝚽lH[:, l] - 𝚽H[:, l], plotOrder = ord, ms = 4)
    # plot!(xlabel = "x(μm)", ylabel = "y(μm)", frame = :box, cbar = true)
    # plot!(left_margin = 3mm, bottom_margin = 3mm, right_margin = 3mm)
    # plt = plot!(size = (500, 500))
    savefig(plt, "../figs/RGC100_HGLET_LPHGLET_j$(j-1)k$(k-1)l$(l-1).png")
end
