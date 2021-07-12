# Script to generate figures in Figure 2.
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms

# First define a utility function "Glevel", which generate a Graph Signal that 
# are piecewise constants whose values are actual region indices "k".
function Glevel(G::GraphSig, GP::GraphPart, j::Int64)
    f = zeros(size(G.f))
    for k in 1:size(GP.rs,1)
        a = GP.rs[k,j+1]
        b = GP.rs[k+1,j+1] - 1
        if b == -1
            break
        end
        f[a:b] .= k*1.0
    end
    Gsub = deepcopy(G)
    Gsub.f[GP.ind] = f
    return Gsub
end

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the Toronto street network (vehicular volume counts as its graph signal)
JLD2.@load "Toronto_new.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],
    plotspecs = tmp1["plotspecs"])

# Assign correct edge weights via 1/Euclidean distance
G = Adj2InvEuc(G)
# Perform the hierarchical graph partitioning using the Fiedler vectors
GP = partition_tree_fiedler(G) # Lrw by default
# Compute the expansion coefficients of the full GHWT c2f dictionary
dmatrix = ghwt_analysis!(G, GP=GP)
N = length(G.f)

# Fig. 2a: Level 1 first since Level 0 is not interesting
j = 1
GraphSig_Plot(Glevel(G,GP,j), linewidth = 1., markersize = 4., 
    markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(showaxis = false, ticks = false)
savefig("toronto_j1.pdf")
savefig("toronto_j1.png")

# Fig. 2b: Level 2
j = 2
GraphSig_Plot(Glevel(G,GP,j), linewidth = 1., markersize = 4., 
    markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(showaxis = false, ticks = false)
savefig("toronto_j2.pdf")
savefig("toronto_j2.png")

# Fig. 2c: Level 3
j = 3
GraphSig_Plot(Glevel(G,GP,j), linewidth = 1., markersize = 4., 
    markercolor = :viridis, markerstrokealpha =0., notitle = true, nocolorbar = true)
plot!(showaxis = false, ticks = false)
savefig("toronto_j3.pdf")
savefig("toronto_j3.png")
