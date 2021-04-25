cd(@__DIR__); include("setups/grid7x3.jl");

## ROT2
D = natural_eigdist(𝚽, 𝛌, Q; α = 0.1, distance = :ROT)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
E[1, :] .*= sign(E[1, 2]); E[2, :] .*= -sign(E[1, 2])
plt = grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot)
savefig(plt, "../figs/Grid7x3_MDS_ROT2_alpha01.png")
display(plt)
