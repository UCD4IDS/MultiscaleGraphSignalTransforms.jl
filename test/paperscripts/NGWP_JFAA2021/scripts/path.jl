# script for Fig.5

using MultiscaleGraphSignalTransforms, LightGraphs, Plots; gr(dpi = 200)
import WaveletsExt: wiggle

## Build Graph
N = 512; G = path_graph(N)
X = zeros(N,2); X[:, 1] = 1:N
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1,:])'
W = 1.0 * adjacency_matrix(G)

## Build NGWPs
Gstar_Sig = GraphSig(W)
G_Sig = GraphSig(W, xy = X)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

@time VM_NGWP = vm_ngwp(𝚽, GP_dual)

#################### Fig.5
j = 5
for k in [1, 2, 5]
   WW = sort_wavelets(VM_NGWP[GP_dual.rs[k, j]:(GP_dual.rs[k + 1, j] - 1), j, :]')
   if k == 2
      WW[:, end] *= -1
   end
   plt = wiggle(WW; sc = 0.75)
   # savefig(plt, joinpath(@__DIR__, "../paperfigs/Path512_VM_NGWP_j$(j-1)k$(k-1).png"))
end
current()
