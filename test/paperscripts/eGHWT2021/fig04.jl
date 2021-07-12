# Script to generate figures in Figure 4.
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms

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

# Fig. 4a: Sacling vector ψ^1_{0,0}, i.e., (j,k,l)=(1,0,0).
j = 1
loc = 1
BS = bs_level(GP, j)
dvec = zeros(N, 1)
dvec[loc,1] = 1.0
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha = 0., 
    markercolor = :viridis, notitle = true, nocolorbar = true, clim = (-0.01, 0.01))
plot!(showaxis = false, ticks = false)
k, l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
savefig("toronto_psi100.pdf")
savefig("toronto_psi100.png")

# Fig. 4b: Haar vector ψ^2_{0,1}, i.e., (j,k,l)=(2,0,1).
j = 2
loc = 2
BS = bs_level(GP, j)
dvec = zeros(N, 1)
dvec[loc,1] = 1.0
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha = 0., 
    markercolor = :viridis, notitle = true, nocolorbar = true, clim = (-0.01, 0.01))
plot!(showaxis = false, ticks = false)
k, l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
savefig("toronto_psi201.pdf")
savefig("toronto_psi201.png")

# Fig. 4c: Walsh vector ψ^3_{0,9}, i.e., (j,k,l)=(3,0,9).
#j = 5
#loc = 1000
j = 3
loc = 10
BS = bs_level(GP, j)
dvec = zeros(N, 1)
dvec[loc,1] = 1
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha = 0., 
    markercolor = :viridis, notitle = true, nocolorbar = true, clim = (-0.01, 0.01))
plot!(showaxis = false, ticks = false)
k, l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
savefig("toronto_psi309.pdf")
savefig("toronto_psi309.png")

