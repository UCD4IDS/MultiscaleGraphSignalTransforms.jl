using NGWP, LightGraphs, MTSG, Plots, LaTeXStrings

## Build Graph
N = 512; G = path_graph(N)
X = zeros(N,2); X[:, 1] = 1:N
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1,:])'
W = 1.0 * adjacency_matrix(G)

Gstar_Sig = GraphSig(W)
G_Sig = GraphSig(W, xy = X)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

## utility functions
function find_mainsupport(w; ϵ = 0.01)
   N = length(w)
   l, r = 1, N
   for i in 1:N
      if abs(w[i]) >= ϵ
         l = i
         break
      end
   end
   for i in N:-1:l
      if abs(w[i]) >= ϵ
         r = i
         break
      end
   end
   return [l, r]
end

function findlocalmaxima(signal::Vector)
    inds = Int[]
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
               push!(inds,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end
    inds
end

function sidelobe_attenuation(w)
   locmax_ind = findlocalmaxima(w)
   locmax_val = sort(w[locmax_ind]; rev = true)
   return locmax_val[2] / locmax_val[1]
end
