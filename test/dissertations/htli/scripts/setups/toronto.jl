using MultiscaleGraphSignalTransforms, JLD, Plots, Graphs, Distances

## Build weighted toronto street network graph
G = loadgraph("../datasets/new_toronto_graph.lgz"); N = nv(G)
X = load("../datasets/new_toronto.jld", "xy")
dist_X = pairwise(Euclidean(), X; dims = 1)
A = 1.0 .* adjacency_matrix(G)
W = zeros(N, N); W[A .> 0] = 1 ./ dist_X[A .> 0]; W = A .* W
Q = incidence_matrix(G; oriented = true)
edge_weight = 1 ./ sqrt.(sum((Q' * X).^2, dims = 2)[:])

## eigenvectors of L(G)
deg = sum(W, dims = 1)[:]  # weighted degree vector
L = diagm(deg) - W
𝛌, 𝚽 = eigen(L)
standardize_eigenvectors!(𝚽)

## eigenvectors of Lsym(G)
Lsym = diagm(deg.^(-1/2)) * (diagm(deg) - W) * diagm(deg.^(-1/2))
𝛌sym, 𝚽sym = eigen(Lsym)
standardize_eigenvectors!(𝚽sym)

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽, Q, N; edge_weight = edge_weight) #52.375477 seconds
Gstar_Sig = dualgraph(distDAG)
G_Sig = GraphSig(A, xy = X); G_Sig = Adj2InvEuc(G_Sig)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)
jmax = size(GP_dual.rs, 2) - 1  # zero-indexed

if allNGWPs
    #1315.821724 seconds (3.05 M allocations: 495.010 GiB, 7.04% gc time)
    @time VM_NGWP = vm_ngwp(𝚽, GP_dual)
    #119.590168 seconds (12.14 M allocations: 158.035 GiB, 13.89% gc time)
    @time PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal)
    @time LP_NGWP = lp_ngwp(𝚽, Gstar_Sig.W, GP_dual; ϵ = 0.3)
end


## Build Dual Graph by DAG metric (Lsym)
distDAG_Lsym = eigDAG_Distance(𝚽sym, Q, N; edge_weight = edge_weight)
Gstar_Sig_Lsym = dualgraph(distDAG_Lsym)
GP_dual_Lsym = partition_tree_fiedler(Gstar_Sig_Lsym; swapRegion = false)
GP_primal_Lsym = pairclustering(𝚽sym, GP_dual_Lsym)
jmax_Lsym = size(GP_dual_Lsym.rs, 2) - 1

if allNGWPs
    VM_NGWP_Lsym = vm_ngwp(𝚽sym, GP_dual_Lsym)
    PC_NGWP_Lsym = pc_ngwp(𝚽sym, GP_dual_Lsym, GP_primal_Lsym)
    LP_NGWP_Lsym = lp_ngwp(𝚽sym, Gstar_Sig_Lsym.W, GP_dual_Lsym; ϵ = 0.3)
end
