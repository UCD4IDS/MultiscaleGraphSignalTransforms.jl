cd(@__DIR__); allNGWPs = false; include("setups/toronto.jl")
gr(dpi = 200)

## pedestrian signal 16 most important LP-NGWP vectors
LP_NGWP = lp_ngwp(𝚽, Gstar_Sig.W, GP_dual; ϵ = 0.3)
fp = load("../datasets/new_toronto.jld", "fp")
G_Sig.f = reshape(fp, (N, 1))
dmatrix_LP = ngwp_analysis(G_Sig, LP_NGWP)
dvec_lp_ngwp, BS_lp_ngwp = ngwp_bestbasis(dmatrix_LP, GP_dual)
important_idx = sortperm(dvec_lp_ngwp[:].^2; rev = true)
println("================fp-LP-NGWP-top-basis-vectors=================")
for i in 1:16
    dr, dc = BS_lp_ngwp.levlist[important_idx[i]]
    w = LP_NGWP[dr, dc, :]
    j, k, l = NGWP_jkl(GP_dual, dr, dc)
    println("(j, k, l) = ($(j), $(k), $(l))")
    sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(A, X; width=1)
    scatter_gplot!(X; marker = sgn .* w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    savefig(plt,
        "../figs/Toronto_fp_DAG_LP_NGWP_ibv$(lpad(i,2,"0"))_j$(j)_k$(k)_l$(l).png")
end
