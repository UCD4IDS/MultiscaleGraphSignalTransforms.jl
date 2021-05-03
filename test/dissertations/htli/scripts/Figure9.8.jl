cd(@__DIR__); allNGWPs = false; include("setups/toronto.jl")
using Plots.PlotMeasures
gr(dpi = 200)

## (a) pedestrian volume graph signal
fp = load("../datasets/new_toronto.jld", "fp")
G_Sig.f = reshape(fp, (N, 1))
gplot(A, X; width=1); plot!(size = (600, 600), right_margin = 5mm)
signal_plt = scatter_gplot!(X; marker = fp, plotOrder = :s2l, ms = 3)
savefig(signal_plt, "../figs/Toronto_fp.png")

## (b) pedestrian signal relative l2 approximation error by various methods
# uncomment the following if `allNGWPs` = true
# DVEC = getall_expansioncoeffs2(G_Sig, GP_dual, GP_dual_Lsym, VM_NGWP, PC_NGWP, LP_NGWP,
#                                VM_NGWP_Lsym, PC_NGWP_Lsym, LP_NGWP_Lsym, 𝚽, 𝚽sym)
# use precomputed results
DVEC = load("../datasets/Toronto_fp_DVEC.jld", "DVEC")
approx_error_plot2(DVEC);
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 11, size = (600, 600))
savefig(plt, "../figs/Toronto_fp_DAG_approx.png")
