using Plots, Graphs, JLD, LaTeXStrings, MultiscaleGraphSignalTransforms, Distances
using Plots.PlotMeasures

G = loadgraph("../datasets/RGC100.lgz"); N = nv(G)
X = load("../datasets/RGC100_xyz.jld", "xyz")[:, 1:2]
X3 = load("../datasets/RGC100_xyz.jld", "xyz")
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L)
standardize_eigenvectors!(𝚽)

dist_X = pairwise(Euclidean(1e-12), X3; dims = 1)
A = 1.0 .* adjacency_matrix(G)
W = zeros(N, N); W[A .> 0] = 1 ./ dist_X[A .> 0]; W = A .* W

Q = incidence_matrix(G; oriented = true)
∇𝚽 = Q' * 𝚽
edge_length = sqrt.(sum((Q' * X3).^2, dims = 2)[:])

G_Sig = GraphSig(W)
GP = partition_tree_fiedler(G_Sig; swapRegion = false)
