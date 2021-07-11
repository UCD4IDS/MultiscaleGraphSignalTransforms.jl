# Script to generate Figures 16 and 17
using TestImages, Plots, SparseArrays, LinearAlgebra, Wavelets, MultiscaleGraphSignalTransforms
include("auxilaries.jl")

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the cameraman image and subsample it to make it as 128x128 image
img = testimage("camera");
matrix = convert(Array{Float64,2}, img)[1:4:512,1:4:512]
m, n = size(matrix)

# Fig. 16a
display(heatmap(matrix, ratio=1, yaxis =:flip, showaxis = false, ticks = false,
    color = :grays, clim = (0,1), colorbar = false))
savefig("cameraman.pdf")
savefig("cameraman.png")

# Now, let's compute the weight matrix based on the Gaussian affinity
# of the small windows (or radius r) around pixels and the pixel locations.
# Set up the key parameters
r = 5; # specify radius of neighbors
σ = 0.007
# Do the weight matrix computation
W007 = image_Gaussian_affinity(matrix, r, σ)    

# Preprocess to generate G (GraphSig struct) and GP (GraphPart struct)
# Note that GraphSig requires a matrix data even if it is just a one vector, i.e.,
# f = matrix[:] does not work!
G007 = GraphSig(W007, f = reshape(matrix, (length(matrix), 1)))
GP007 = partition_tree_fiedler(G007)
dmatrix007 = ghwt_analysis!(G007, GP=GP007)

# Construct or search the specific basis
# Haar
BS_haar007 = bs_haar(GP007)
dvec_haar007 = dmatrix2dvec(dmatrix007, GP007, BS_haar007)

# eGHWT
dvec_eghwt007, BS_eghwt007 = eghwt_bestbasis(dmatrix007, GP007)

# Generate Fig. 17a
top_vectors_plot2(dvec_eghwt007, m, n, BS_eghwt007, GP007)
savefig("cameraman_eghwt09_sigma007.pdf")
savefig("cameraman_eghwt09_sigma007.png")

# Now change the σ value
σ = 0.07
W07 = image_Gaussian_affinity(matrix, r, σ)    
# Preprocess to generate G (GraphSig struct) and GP (GraphPart struct)
G07 = GraphSig(W07, f = reshape(matrix, (length(matrix), 1)))
GP07 = partition_tree_fiedler(G07)
dmatrix07 = ghwt_analysis!(G07, GP=GP07)

# Construct or search the specific basis
# Haar
BS_haar07 = bs_haar(GP07)
dvec_haar07 = dmatrix2dvec(dmatrix07, GP07, BS_haar07)
# eGHWT
dvec_eghwt07, BS_eghwt07 = eghwt_bestbasis(dmatrix07, GP07)
# Generate Fig. 17b
top_vectors_plot2(dvec_eghwt07, m, n, BS_eghwt07, GP07)
savefig("cameraman_eghwt09_sigma07.pdf")
savefig("cameraman_eghwt09_sigma07.png")
# Finally, do the classical Haar transform in Wavelets.jl
dvec_classichaar = dwt(matrix, wavelet(WT.haar))

# Fig. 16b (Approximation error plot)
DVEC = [ dvec_classichaar[:], dvec_haar07[:], dvec_eghwt07[:],
         dvec_haar007[:], dvec_eghwt007[:] ]
T = [ "Classical Haar", "Graph Haar (σ = 0.07)", 
      "eGHWT (σ = 0.07)", "Graph Haar (σ = 0.007)", 
      "eGHWT (σ = 0.007)"]
L = [ (:dashdot,:orange), (:dashdot, :red), (:dashdot, :black),
      (:solid, :red), (:solid, :black) ]
approx_error2(DVEC, T, L, 0.5)
savefig("cameraman_approx_error.pdf")
savefig("cameraman_approx_error.png")
