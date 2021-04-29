cd(@__DIR__); include("setups/grid7x3.jl");

## ROT1 + pmf1
D = natural_eigdist(𝚽, 𝛌, Q; α = 0.5, input_format = :pmf1, distance = :ROT)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
plt = grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot)
savefig(plt, "../figs/Grid7x3_MDS_ROT1_pmf1_alpha05.png")
display(plt)
