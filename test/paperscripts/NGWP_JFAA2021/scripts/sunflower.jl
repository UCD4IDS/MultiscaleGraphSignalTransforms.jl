# script for Fig.8(b)(c), Fig.9, Fig.10(b)(c), Fig.11

using MultiscaleGraphSignalTransforms, LightGraphs, Plots; gr(dpi = 200)

## Build weighted graph
G, L, X = SunFlowerGraph(N = 400); N = nv(G)
𝛌, 𝚽 = eigen(Matrix(L))
sgn = (maximum(𝚽, dims = 1)[:] .> -minimum(𝚽, dims = 1)[:]) .* 2 .- 1
𝚽 = 𝚽 * Diagonal(sgn)
Q = incidence_matrix(G; oriented = true)
W = 1.0 * adjacency_matrix(G)
edge_weight = [e.weight for e in edges(G)]

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽, Q, N; edge_weight = edge_weight)
Gstar_Sig = dualgraph(distDAG)
G_Sig = GraphSig(W, xy = X)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

@time VM_NGWP = vm_ngwp(𝚽, GP_dual) #54.986524 seconds (1.19 M allocations: 25.127 GiB, 2.36% gc time)
@time PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal) #0.611176 seconds (225.41 k allocations: 844.488 MiB, 11.71% gc time)


#################### Fig. 8(b) barbara eye graph signal
using MAT
f = matread(joinpath(@__DIR__, "../datasets",
            "sunflower_barbara_voronoi.mat"))["f_eye_voronoi"]
G_Sig.f = reshape(f, (N, 1))
scatter_gplot(X; marker = f, ms = LinRange(4.0, 14.0, N), c = :greys);
plt = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], frame = :none)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/SunFlower_barbara_feye.png"))

#################### Fig. 8(c) barbara eye relative l2 approximation error by various methods
DVEC = getall_expansioncoeffs(G_Sig, GP_dual, VM_NGWP, PC_NGWP, 𝚽)
approx_error_plot(DVEC);
plt = plot!(xguidefontsize = 16, yguidefontsize = 16, legendfontsize = 12)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/SunFlower_barbara_feye_DAG_approx.png"))

#################### Fig. 9 barbara eye 16 most important VM-NGWP vectors (ignore the DC vector)
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
for i in 2:17
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    println("(j, k, l) = ", NGWP_jkl(GP_dual, dr, dc))
    scatter_gplot(X; marker = w, ms = LinRange(4.0, 14.0, N), c = :greys)
    plt = plot!(xlim = [-1.2, 1.2], ylim = [-1.2, 1.2], frame = :none,
                cbar = false, clims = (-0.15, 0.15))
    # savefig(plt, joinpath(@__DIR__,
    #     "../paperfigs/SunFlower_barbara_feye_DAG_VM_NGW_important_basis_vector$(lpad(i,2,"0")).png"))
end

#################### Fig. 10(b) barbara pants graph signal
f = matread(joinpath(@__DIR__, "../datasets",
            "sunflower_barbara_voronoi.mat"))["f_trouser_voronoi"]
scatter_gplot(X; marker = f, ms = LinRange(4.0, 14.0, N), c = :greys);
plt = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], frame = :none)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/SunFlower_barbara_ftrouser.png"))

#################### Fig. 10(c) barbara eye relative l2 approximation error by various methods
G_Sig.f = reshape(f, (N, 1))
DVEC = getall_expansioncoeffs(G_Sig, GP_dual, VM_NGWP, PC_NGWP, 𝚽)
approx_error_plot(DVEC);
plt = plot!(xguidefontsize = 16, yguidefontsize = 16, legendfontsize = 12)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/SunFlower_barbara_ftrouser_DAG_approx.png"))

#################### Fig. 11 barbara pants 16 most important VM-NGWP vectors (ignore the DC vector)
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
for i in 2:17
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    println("(j, k, l) = ", NGWP_jkl(GP_dual, dr, dc))
    scatter_gplot(X; marker = w, ms = LinRange(4.0, 14.0, N), c = :greys)
    plt = plot!(xlim = [-1.2, 1.2], ylim = [-1.2, 1.2], frame = :none,
                cbar = false, clims = (-0.15, 0.15))
    # savefig(plt, joinpath(@__DIR__,
    #     "../paperfigs/SunFlower_barbara_ftrouser_DAG_VM_NGW_important_basis_vector$(lpad(i,2,"0")).png"))
end
