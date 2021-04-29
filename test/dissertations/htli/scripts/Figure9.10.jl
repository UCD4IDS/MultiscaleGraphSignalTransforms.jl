cd(@__DIR__); allNGWPs = false; include("setups/toronto.jl")
gr(dpi = 200)

## pedestrian signal 16 most important PC-NGWP vectors
PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal)
fp = load("../datasets/new_toronto.jld", "fp")
G_Sig.f = reshape(fp, (N, 1))
dmatrix_PC = ngwp_analysis(G_Sig, PC_NGWP)
dvec_pc_ngwp, BS_pc_ngwp = ngwp_bestbasis(dmatrix_PC, GP_dual)
important_idx = sortperm(dvec_pc_ngwp[:].^2; rev = true)
println("================fp-PC-NGWP-top-basis-vectors=================")
for i in 1:16
    dr, dc = BS_pc_ngwp.levlist[important_idx[i]]
    w = PC_NGWP[dr, dc, :]
    j, k, l = NGWP_jkl(GP_dual, dr, dc)
    print("(j, k, l) = ($(j), $(k), $(l))  ")
    if j == jmax
        print("φ_{$(GP_dual.ind[dr]-1)}")
    end
    println()
    sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(A, X; width=1)
    scatter_gplot!(X; marker = sgn .* w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    savefig(plt,
        "../figs/Toronto_fp_DAG_PC_NGWP_ibv$(lpad(i,2,"0"))_j$(j)_k$(k)_l$(l).png")
end
