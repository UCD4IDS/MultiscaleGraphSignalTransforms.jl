# Script to generate Figure 15.
using Plots, SparseArrays, JLD2, LinearAlgebra, Wavelets, MultiscaleGraphSignalTransforms
include("auxilaries.jl")
include("../../../src/unbalanced_haar_image.jl")

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the Barbara image data
JLD2.@load "barbara.jld2"

# smaller face region of size 100x100
row_zoom = 71:170
col_zoom = 341:440
display(heatmap(barbara[row_zoom, col_zoom],ratio=1, yaxis =:flip, 
    showaxis = false, ticks = false, color = :grays, clims = (0,1)))

# Extract that portion of the image
matrix = deepcopy(barbara[row_zoom, col_zoom])

#
# Perform the image partition using the penalized total variation and then
# compute the expansion coefficients
#
p = 3.0
maxlev = 9 # need this to reach the single node leaves
GPcols = PartitionTreeMatrix_unbalanced_haar_binarytree(matrix, maxlev, p)
maxlev = 8 # need this to reach the single node leaves
GProws = PartitionTreeMatrix_unbalanced_haar_binarytree(copy(matrix'), maxlev, p)
# copy() is necessary to switch the Adjoint type to a regular Matrix{Float64}.
dmatrix = ghwt_analysis_2d(matrix, GProws, GPcols)

#
# Compute the graph Haar coefficients from the previous PTV partition tree
#
BS_haar_rows = bs_haar(GProws)
BS_haar_cols = bs_haar(GPcols)
dvec_haar, loc_haar = BS2loc(dmatrix, GProws, GPcols, BS_haar_rows, BS_haar_cols)

#
# Compute the eGHWT best basis
#
dvec_eghwt, loc_eghwt = eghwt_bestbasis_2d(matrix, GProws, GPcols)

#
# Classical Haar transform via Wavelets.jl with direct input
#
dvec_classichaar = dwt(matrix, wavelet(WT.haar))

#
# Classical Haar transform via Wavelets.jl via putting in the dyadic image
#
matrix2 = zeros(128,128)
matrix2[15:114,15:114] = deepcopy(matrix)
dvec_classichaar2 = dwt(matrix2, wavelet(WT.haar))

#
# Classical Haar transform via Wavelets.jl via putting in the dyadic image
#
matrix3 = zeros(128,128)
matrix3[1:100,1:100] = deepcopy(matrix)
matrix3[101:end,1:100] = deepcopy(matrix[end:-1:end-27,1:100])
matrix3[:,101:end] = matrix3[:,100:-1:73]
dvec_classichaar3 = dwt(matrix3, wavelet(WT.haar))

# Figure 15 (Approximation error plot up to the top 5000 coefficients)
DVEC = [ dvec_classichaar[:], dvec_classichaar2[:], dvec_classichaar3[:], 
    dvec_haar[:], dvec_eghwt[:] ]
T = [ "Classic Haar (direct)", "Classic Haar (zero pad)", 
    "Classic Haar (even refl)", "Graph Haar (PTV cost)", "eGHWT (PTV cost)" ]
L = [ (:dashdot, :orange), (:dashdot, :red), (:dashdot, :green),
      (:solid, :blue), (:solid, :black) ]
approx_error3(DVEC, T, L, 5000)
savefig("barbara_face100_nondyadic_haar_error.pdf")
savefig("barbara_face100_nondyadic_haar_error.png")
