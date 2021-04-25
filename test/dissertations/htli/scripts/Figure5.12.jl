cd(@__DIR__); include("setups/simpletree.jl");
plotlyjs(dpi = 200)

## (a)
D = natural_eigdist(𝚽, 𝛌, Q; α = 1.0, input_format = :pmf1, distance = :ROT)
E = transform(fit(MDS, D, maxoutdim=3, distances=true))
E[1, :] .*= 1; E[2, :] .*= -1; E[3, :] .*= 1;
plt = simpletree_mds_plot(E, ieb1, ieb2, ieb3, ieb4, iejc); plot!(camera=(210,-37), xlims = [-18, 15])
# savefig(plt, "../figs/SimpleTree_MDS_ROT1_pmf1_alpha1.html")
savefig(plt, "../figs/SimpleTree_MDS_ROT1_pmf1_alpha1.pdf")

## (b)
D = dist_sROT
E = transform(fit(MDS, D, maxoutdim=3, distances=true))
plt = simpletree_mds_plot(E, ieb1, ieb2, ieb3, ieb4, iejc)
# savefig(plt, "../figs/SimpleTree_MDS_sROT_pmf1_alpha1.html")
savefig(plt, "../figs/SimpleTree_MDS_sROT_pmf1_alpha1.pdf")
