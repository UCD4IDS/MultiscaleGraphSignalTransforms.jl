using MultiscaleGraphSignalTransforms, Graphs, Plots, LaTeXStrings

## Build Graph
N = 16; G = path_graph(N)
X = zeros(N,2); X[:, 1] = 1:N
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1,:])'
W = 1.0 * adjacency_matrix(G)

## preliminaries for PC-NGWP
Gstar_Sig = GraphSig(W)
G_Sig = GraphSig(W, xy = X)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

## Plots
function path_pc_plot!(W, X; ind = 1:size(X, 1), c = :teal)
    # plot(size = (500, 50), framestyle = :none, xlim = [1, N+1.5], ylim = [-0.2, 0.3])
    gplot!(W, X, width = 1, color = :blue, style = :solid)
    scatter_gplot!(X; ms = 6, c = :teal)
    annotate!(X[:, 1], X[:, 2] .- 0.2, [text(string(n), :top, 7) for n = 1:N])
    scatter_gplot!(X[ind, :]; ms = 6, c = c)
    plt = plot!(ratio = :auto)
    return plt
end
