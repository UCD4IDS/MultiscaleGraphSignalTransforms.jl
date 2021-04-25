cd(@__DIR__); include("setups/grid7x3.jl");

## ROT2
D = natural_eigdist(𝚽, 𝛌, Q; α = 1.0, distance = :ROT)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
plt = grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot,
            annotate_ind = vcat(1:6, 8, 10, 11), plotOrder = N:-1:1)
savefig(plt, "../figs/Grid7x3_MDS_ROT2_alpha1.png")
display(plt)

## TSD (T = :Inf)
D = natural_eigdist(𝚽, 𝛌, Q; T = :Inf, distance = :TSD)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
E[1, :] .*= sign(E[1, 2]); E[2, :] .*= -sign(E[1, 2])
plt = grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot,
            annotate_ind = vcat(1:8, 11), plotOrder = vcat(N:-1:9, 6:8, 5:-1:1))
savefig(plt, "../figs/Grid7x3_MDS_TSD_Tinfty.png")
display(plt)


## TSD T from 0.1 to 10
anim = @animate for t ∈ 0.1:0.1:10
    D = natural_eigdist(𝚽, 𝛌, Q; T = t, distance = :TSD)
    E = transform(fit(MDS, D, maxoutdim=2, distances=true))
    E[1, :] .*= sign(E[1, 2])
    E[2, :] .*= sign(E[2, 2]) * ((t <= 1) * 2 - 1)
    grid7x3_mds_heatmaps(E, 𝚽; backend = :pyplot)
    title!("T = $(t)")
end
gif(anim, "../gifs/Grid7x3_MDS_TSD_T01_inf.gif", fps = 5)
