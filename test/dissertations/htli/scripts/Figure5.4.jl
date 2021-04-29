cd(@__DIR__); include("setups/grid7x3.jl");

## HAD
D = natural_eigdist(𝚽, 𝛌, Q; distance = :HAD)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
plt = grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot)
savefig(plt, "../figs/Grid7x3_MDS_HAD.png")
display(plt)
