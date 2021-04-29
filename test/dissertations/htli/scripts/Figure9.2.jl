cd(@__DIR__); runapprox = true; allNGWPs = true; include("setups/sunflower.jl")
gr(dpi = 200)

## (a) barbara sunflower eye
heatmap(barbara, yflip = true, ratio = 1, c = :greys)
scatter_gplot!(transform2D(X; s = 20, t = [395, 100]); ms = 2, c = :red)
sample_location_plt = plot!(cbar = false, frame = :none)
savefig(sample_location_plt, "../figs/barb_sunflower_eye.png")

## (b) barbara eye graph signal
f = matread("../datasets/sunflower_barbara_voronoi.mat")["f_eye_voronoi"]
G_Sig.f = reshape(f, (N, 1))
scatter_gplot(X; marker = f, ms = LinRange(4.0, 14.0, N), c = :greys);
plt = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], frame = :none)
savefig(plt, "../figs/SunFlower_barbara_feye.png")

## (c) barbara eye relative l2 approximation error by various methods
DVEC = getall_expansioncoeffs2(G_Sig, GP_dual, GP_dual_Lsym, VM_NGWP, PC_NGWP, LP_NGWP,
                               VM_NGWP_Lsym, PC_NGWP_Lsym, LP_NGWP_Lsym, 𝚽, 𝚽sym)
approx_error_plot2(DVEC);
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 11)
savefig(plt, "../figs/SunFlower_barbara_feye_DAG_approx.png")
