# Generate figures on the Toronto street network

# Load necessary packages
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms
include("auxilaries.jl")

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the Toronto street network (vehicular volume counts as its graph signal)
JLD2.@load "Toronto_new.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],
    plotspecs = tmp1["plotspecs"])

# Generate Fig. 10a
# GraphSig_Plot(G, linewidth = 1., markersize = 4., 
#    markercolor = :viridis, markerstrokealpha = 0., notitle = true)
gplot(G.W, G.xy; width=1); 
plt = scatter_gplot!(G.xy; marker = G.f, plotOrder = :s2l, 
    ms = G.f*8/maximum(G.f), mswidth=1)
display(plt)
savefig("toronto_vv.pdf")
savefig("toronto_vv.png")

# Assign correct edge weights via 1/Euclidean distance
G = Adj2InvEuc(G)
# Perform the hierarchical graph partitioning using the Fiedler vectors
GP = partition_tree_fiedler(G) # Lrw by default
# Compute the expansion coefficients of the full GHWT c2f dictionary
dmatrix = ghwt_analysis!(G, GP=GP)


# Extract the coefficients corresponding to the Haar basis
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

# Extract the coefficients corresponding to the Walsh basis
BS_walsh = bs_walsh(GP)
dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)

# Compute the c2f GHWT best basis and the coefficients
dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)

# Compute the f2c GHWT best basis and the coefficients
dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)

# Compute the eGHWT best basis and the coefficients
dvec_eghwt, BS_eghwt = eghwt_bestbasis(dmatrix, GP)

################################################################################
###################### Generate approximaation errors ##########################
################################################################################
# Generate Fig. 10b up to 50% of the coefficients retained.
DVEC = [ dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:] ]
T = [ "Graph Haar","Graph Walsh","GHWT_c2f", "GHWT_f2c", "eGHWT" ]
L = [ (:dashdot, :orange), (:dashdot, :blue), (:solid, :red), (:solid, :green),
      (:solid, :black) ]
approx_error2(DVEC, T, L, 0.5)
savefig("toronto_vv_approx_error.pdf")
savefig("toronto_vv_approx_error.png")


################################################################################
######################### Synthesis by top vectors #############################
################################################################################
color_limit_residual = (0., 0.15)

function top_vectors_residual(p::Int64, dvec::Array{Float64,2}, 
                                BS::BasisSpec, GP::GraphPart, G::GraphSig)
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0
    (recon, GS)= ghwt_synthesis(dvecT, GP, BS, G)
    GS.name = "Square of the residual"
    #GS.f =  (G.f - GS.f).^2
    #print((maximum(GS.f), minimum(GS.f)))
    #GraphSig_Plot(GS)
    GS.f = (G.f - recon).^2 ./(G.f.^2)
    GraphSig_Plot(GS, linewidth = 1., markersize = 4., markercolor = :viridis, 
        markerstrokealpha =0., clim = color_limit_residual)
    display(current())
end

################################################################################
############## Visuallization of Top basis vectors #############################
################################################################################

### function to plot the basis vectors of user specified range
# Recommended number of basis vectors are squared numbers, e.g., 9, 16, 25.
function top_vectors_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart,
                             G::GraphSig; istart = 2, iend = 17, ms = 2)
    sorted_ind = sortperm(dvec[:].^2; rev = true)
    ibs = iend-istart+1
    n1 = Int64(sqrt(ibs))
    plot(layout = Plots.grid(n1, n1))
    idsp = 1
    for ib in istart:iend
        dvecT = fill(0., size(dvec))
        dvecT[sorted_ind[ib]] = 1
        bv = ghwt_synthesis(dvecT, GP, BS)
        gplot!(G.W, G.xy; width=0.25, subplot=idsp, color=:gray)
        scatter_gplot!(G.xy; marker = bv, plotOrder = :propabs, ms, subplot=idsp)
        plt = plot!(cbar = false, clims = (-0.001,0.001), axis=([], false), subplot=idsp)
        idsp += 1
    end
    display(current())
end # of top_vectors_plot

### Fig. 11a: Haar
top_vectors_plot(dvec_haar, BS_haar, GP, G, iend=10)
savefig("toronto_vv_haar09.png")
savefig("toronto_vv_haar09.pdf")

### Walsh
top_vectors_plot(dvec_walsh, BS_walsh, GP, G, iend=10)
savefig("toronto_vv_walsh09.png")
savefig("toronto_vv_walsh09.pdf")

### Fig. 11b: GHWT_c2f (= Walsh in this case)
top_vectors_plot(dvec_c2f, BS_c2f, GP, G, iend=10)
savefig("toronto_vv_c2f09.png")
savefig("toronto_vv_c2f09.pdf")

### Fig. 11c: GHWT_f2c
top_vectors_plot(dvec_f2c, BS_f2c, GP, G, iend=10)
savefig("toronto_vv_f2c09.png")
savefig("toronto_vv_f2c09.pdf")

### Fig. 11d: eGHWT
top_vectors_plot(dvec_eghwt, BS_eghwt, GP, G, iend=10)
savefig("toronto_vv_eghwt09.png")
savefig("toronto_vv_eghwt09.pdf")

