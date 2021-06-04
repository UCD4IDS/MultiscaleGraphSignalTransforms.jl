#Figure 2.5
#Run function "Glevel" at bottom of file first
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms

JLD2.@load "../data/Toronto.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])

G = Adj2InvEuc(G)
GP = partition_tree_fiedler(G,:Lrw)
dmatrix = ghwt_analysis!(G, GP=GP)


j = 1
GraphSig_Plot(Glevel(G,GP,1), linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(axis = false)
#savefig("G1.pdf")

j = 2
GraphSig_Plot(Glevel(G,GP,2), linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(axis = false)
#savefig("G2.pdf")

j = 3
GraphSig_Plot(Glevel(G,GP,3), linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(axis = false)
#savefig("G3.pdf")


function Glevel(G::GraphSig, GP::GraphPart, j::Int64)
    f = zeros(size(G.f))
    for k in 1:size(GP.rs,1)
        a = GP.rs[k,j]
        b = GP.rs[k+1,j] - 1
        if b == -1
            break
        end
        f[a:b] .= k*1.0
    end
    Gsub = deepcopy(G)
    Gsub.f[GP.ind] = f
    return Gsub
end
