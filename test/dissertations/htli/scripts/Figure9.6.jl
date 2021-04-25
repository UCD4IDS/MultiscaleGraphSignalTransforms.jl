cd(@__DIR__); allNGWPs = false; include("setups/toronto.jl")
gr(dpi = 200)

## (a) a smooth spatial distribution of the street intersections graph signal
f = zeros(N); for i in 1:N; f[i] = length(findall(dist_X[:,i] .< 1/minimum(edge_weight))); end
G_Sig.f = reshape(f, (N, 1))
gplot(A, X; width=1); plot!(size = (600, 500))
signal_plt = scatter_gplot!(X; marker = f, plotOrder = :s2l, ms = 3)
savefig(signal_plt, "../figs/Toronto_fdensity.png")

## (b) spatial distribution signal relative l2 approximation error by various methods
# uncomment the following if `allNGWPs` = true
# DVEC = getall_expansioncoeffs2(G_Sig, GP_dual, GP_dual_Lsym, VM_NGWP, PC_NGWP, LP_NGWP,
#                                VM_NGWP_Lsym, PC_NGWP_Lsym, LP_NGWP_Lsym, 𝚽, 𝚽sym)
# use precomputed results
DVEC = load("../datasets/Toronto_fdensity_DVEC.jld", "DVEC")
approx_error_plot2(DVEC);
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 11)
savefig(plt, "../figs/Toronto_fdensity_DAG_approx.png")
