#Figure 3.2
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms
include("../../../../src/utils.jl")

JLD2.@load "../data/Toronto.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])

G = Adj2InvEuc(G)
GP = partition_tree_fiedler(G,:Lrw)
dmatrix = ghwt_analysis!(G, GP=GP)

N = length(G.f)


j = 1
loc = 1
BS = bs_level(GP, j)
dvec = zeros(N,1)
dvec[loc,1] = 1
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha =0., notitle = true, clim = (-0.1, 0.1))
plot!(axis = false,colorbar = false)
k,l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
#savefig("Toronto100.pdf")


j = 2
loc = 2
BS = bs_level(GP, j)
dvec = zeros(N,1)
dvec[loc,1] = 1
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha =0., notitle = true, clim = (-0.1, 0.1))
plot!(axis = false,colorbar = false)
k,l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
#savefig("Toronto201.pdf")


j = 5
loc = 100
BS = bs_level(GP, j)
dvec = zeros(N,1)
dvec[loc,1] = 1
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
GraphSig_Plot(GS, linewidth = 1., markersize = 4., markerstrokealpha =0., notitle = true, clim = (-0.1, 0.1))
plot!(axis = false,colorbar = false)
k,l = BS.levlist[loc][1], BS.levlist[loc][2]
print(j," ",(rs_to_region(GP.rs, GP.tag))[k,l]," ",GP.tag[k,l])
#savefig("Toronto5134.pdf")
