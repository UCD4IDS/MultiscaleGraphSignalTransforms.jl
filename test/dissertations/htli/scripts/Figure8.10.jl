cd(@__DIR__); include("setups/path512.jl")
pyplot(dpi = 200)

## construct full HGLET and Lapped-HGLET dictionaries
@time HGLET_dic = HGLET_dictionary(GP, G_Sig; method = :L)
@time LPHGLET_dic = LPHGLET_dictionary(GP, G_Sig; method = :L, ϵ = 0.3)

j = 3; k = 2; l = 6;
WH = HGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'
WlH = LPHGLET_dic[GP.rs[k, j]:(GP.rs[k + 1, j] - 1), j, :]'

plt = plot(WH[:, l], c = :black, grid = false, frame = :box, lw = 0.8,
     legendfontsize = 11, legend = false, size = (500, 400), ylim = [-0.15, 0.15])
     xticks!([1; 64:64:N], vcat(string(1), [string(k) for k in 64:64:N]))
savefig(plt, "../figs/Path512_HGLET_j$(j-1)k$(k-1)l$(l-1).png")

plt = plot(WlH[:, l], c = :black, grid = false, frame = :box, lw = 0.8,
     legendfontsize = 11, legend = false, size = (500, 400), ylim = [-0.15, 0.15])
     xticks!([1; 64:64:N], vcat(string(1), [string(k) for k in 64:64:N]))
savefig(plt, "../figs/Path512_LPHGLET_j$(j-1)k$(k-1)l$(l-1).png")
