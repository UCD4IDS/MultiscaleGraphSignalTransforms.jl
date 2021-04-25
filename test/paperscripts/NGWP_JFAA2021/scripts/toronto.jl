# script for Fig.12, Fig.13, Fig.14, Fig.15, Fig.16

using MultiscaleGraphSignalTransforms, JLD, Plots, LightGraphs, Distances
gr(dpi = 200)

## Build weighted toronto street network graph
G = loadgraph(joinpath(@__DIR__, "../datasets", "new_toronto_graph.lgz")); N = nv(G)
X = load(joinpath(@__DIR__, "../datasets", "new_toronto.jld"), "xy")
dist_X = pairwise(Euclidean(), X; dims = 1)
W = 1.0 .* adjacency_matrix(G)
Weight = zeros(N, N); Weight[W .> 0] = 1 ./ dist_X[W .> 0]; Weight = W .* Weight
L = Matrix(Diagonal(sum(Weight; dims = 1)[:]) - Weight)
𝛌, 𝚽 = eigen(L);
sgn = (maximum(𝚽, dims = 1)[:] .> -minimum(𝚽, dims = 1)[:]) .* 2 .- 1; 𝚽 = 𝚽 .* sgn';
Q = incidence_matrix(G; oriented = true)
edge_weight = 1 ./ sqrt.(sum((Q' * X).^2, dims = 2)[:])

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽, Q, N; edge_weight = edge_weight) #52.375477 seconds
Gstar_Sig = dualgraph(distDAG)
G_Sig = GraphSig(W, xy = X); G_Sig = Adj2InvEuc(G_Sig)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

@time PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal) #119.590168 seconds (12.14 M allocations: 158.035 GiB, 13.89% gc time)
@time VM_NGWP = vm_ngwp(𝚽, GP_dual) #1315.821724 seconds (3.05 M allocations: 495.010 GiB, 7.04% gc time)

#################### Fig. 12(a) a smooth spatial distribution of the street intersections graph signal
f = zeros(N); for i in 1:N; f[i] = length(findall(dist_X[:,i] .< 1/minimum(edge_weight))); end #fneighbor
G_Sig.f = reshape(f, (N, 1))
gplot(W, X; width=1); signal_plt = scatter_gplot!(X; marker = f, plotOrder = :s2l, ms = 3)
# savefig(signal_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity.png"))

#################### Fig. 12(b) spatial distribution signal relative l2 approximation error by various methods
DVEC = getall_expansioncoeffs(G_Sig, GP_dual, VM_NGWP, PC_NGWP, 𝚽)
approx_error_plot(DVEC);
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 10)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity_DAG_approx.png"))

#################### Fig. 13 fdensity 16 most important VM-NGWP vectors (ignore the DC vector)
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
for i in 2:17
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    println("(j, k, l) = ", NGWP_jkl(GP_dual, dr, dc))
    gplot(W, X; width=1)
    scatter_gplot!(X; marker = w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    # savefig(plt, joinpath(@__DIR__,
    #     "../paperfigs/Toronto_fdensity_DAG_VM_NGW_important_basis_vector$(lpad(i,2,"0")).png"))
end

#################### Fig. 14(a) pedestrian volume graph signal
fp = load(joinpath(@__DIR__, "../datasets", "new_toronto.jld"), "fp")
G_Sig.f = reshape(fp, (N, 1))
gplot(W, X; width=1); signal_plt = scatter_gplot!(X; marker = fp, plotOrder = :s2l, ms = 3)
# savefig(signal_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fp.png"))

#################### Fig. 14(b) pedestrian signal relative l2 approximation error by various methods
DVEC = getall_expansioncoeffs(G_Sig, GP_dual, VM_NGWP, PC_NGWP, 𝚽)
approx_error_plot(DVEC)
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 10)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/Toronto_fp_DAG_approx.png"))

#################### Fig. 15 pedestrian signal 16 most important VM-NGWP vectors
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
for i in 1:16
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    println("(j, k, l) = ", NGWP_jkl(GP_dual, dr, dc))
    gplot(W, X; width=1)
    scatter_gplot!(X; marker = w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    # savefig(plt, joinpath(@__DIR__,
    #     "../paperfigs/Toronto_fp_DAG_VM_NGW_important_basis_vector$(lpad(i,2,"0")).png"))
end

#################### Fig. 16 pedestrian signal 16 most important PC-NGWP vectors
dmatrix_PC = ngwp_analysis(G_Sig, PC_NGWP)
dvec_pc_ngwp, BS_pc_ngwp = ngwp_bestbasis(dmatrix_PC, GP_dual)
important_idx = sortperm(dvec_pc_ngwp[:].^2; rev = true)
for i in 1:16
    dr, dc = BS_pc_ngwp.levlist[important_idx[i]]
    w = PC_NGWP[dr, dc, :]
    println("(j, k, l) = ", NGWP_jkl(GP_dual, dr, dc))
    sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(W, X; width=1)
    scatter_gplot!(X; marker = sgn .* w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    # savefig(plt, joinpath(@__DIR__,
    #     "../paperfigs/Toronto_fp_DAG_PC_NGW_important_basis_vector$(lpad(i,2,"0")).png"))
end
