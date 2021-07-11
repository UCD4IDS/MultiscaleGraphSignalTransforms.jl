# Script to generate Figures 12, 13, and 14.
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms
include("auxilaries.jl")

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the Barbara image data
JLD2.@load "barbara.jld2"
matrix = deepcopy(barbara)

# Here are the list of local regions of interest:
# full image
row_zoom = 1:512
col_zoom = 1:512

# Generate the original Barbara image as a reference.
heatmap(matrix[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, color = :grays)
savefig("barbara.pdf")
savefig("barbara.png")

#
# Initialize the regular balanced binary partition and compute the expansion coefficients
# 
dmatrix, GProws, GPcols = eghwt_init_2d_Lindberg(matrix)

#
# Compute the coefficients of different bases with which we compare
#
# Haar
BS_haar_rows = bs_haar(GProws)
BS_haar_cols = bs_haar(GPcols)
dvec_haar, loc_haar = BS2loc(dmatrix, GProws, GPcols, BS_haar_rows, BS_haar_cols)

# Walsh
BS_walsh_rows = bs_walsh(GProws)
BS_walsh_cols = bs_walsh(GPcols)
dvec_walsh, loc_walsh = BS2loc(dmatrix, GProws, GPcols, BS_walsh_rows, BS_walsh_cols)

# GHWT c2f
fcols, jmax_col = size(GProws.tag);
frows, jmax_row = size(GPcols.tag);
dmatrix_rows = reshape(dmatrix, (frows, jmax_row, fcols*jmax_col))
dmatrix_cols = Array{Float64,3}(reshape(dmatrix',(fcols, jmax_col, frows*jmax_row)))
dvec_c2f_rows, BS_c2f_rows = ghwt_c2f_bestbasis(dmatrix_rows, GProws)
dvec_c2f_cols, BS_c2f_cols = ghwt_c2f_bestbasis(dmatrix_cols, GPcols)
dvec_c2f, loc_c2f = BS2loc(dmatrix, GProws, GPcols, BS_c2f_rows, BS_c2f_cols)

# GHWT f2c
dvec_f2c_rows, BS_f2c_rows = ghwt_f2c_bestbasis(dmatrix_rows, GProws)
dvec_f2c_cols, BS_f2c_cols = ghwt_f2c_bestbasis(dmatrix_cols, GPcols)
dvec_f2c, loc_f2c = BS2loc(dmatrix, GProws, GPcols, BS_f2c_rows, BS_f2c_cols)


# eGHWT
dvec_eghwt, loc_eghwt = eghwt_bestbasis_2d(matrix, GProws, GPcols)

################################################################################
##################### Visualize the results of synthesis #######################
################################################################################

### function to plot the image synthesized by top p vectors
function top_vectors_synthesis_2d(p::Int64, dvec::Vector{Float64}, loc::Matrix{Int64}, 
    GProws::GraphPart, GPcols::GraphPart, dmatrix::Matrix{Float64})
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0

    matrix_syn = eghwt_synthesis_2d(dvecT, loc, GProws, GPcols)
    display(heatmap(matrix_syn[row_zoom, col_zoom],ratio=1, yaxis =:flip, 
        showaxis = false, ticks = false, colorbar= false, 
        clims = (minimum(matrix), maximum(matrix)), color = :grays))
    mse = norm(matrix - matrix_syn,2)^2/length(matrix)
    snr = 20 * log10(norm(matrix,2)/norm(matrix - matrix_syn,2))
    psnr = 10 * log10(maximum(matrix)^2/mse)
    println("MSE: ", mse, "SNR (dB): ", snr, "PSNR (dB): ", psnr)

    return matrix_syn
end

################################################################################
percent = 1/32;
p = Int64(floor(percent*length(matrix)))

# Fig. 13a: Haar
haar32 = top_vectors_synthesis_2d(p, dvec_haar, loc_haar, GProws, GPcols, dmatrix)
savefig("barbara_haar32.pdf")
savefig("barbara_haar32.png")

# Walsh
walsh32 = top_vectors_synthesis_2d(p, dvec_walsh, loc_walsh, GProws, GPcols, dmatrix)
savefig("barbara_walsh32.pdf")
savefig("barbara_walsh32.png")

# Fig. 13b: GHWT c2f
c2f32 = top_vectors_synthesis_2d(p, dvec_c2f, loc_c2f, GProws, GPcols, dmatrix)
savefig("barbara_c2f32.pdf")
savefig("barbara_c2f32.png")

# Fig. 13c: GHWT f2c
f2c32 = top_vectors_synthesis_2d(p, dvec_f2c, loc_f2c, GProws, GPcols, dmatrix)
savefig("barbara_f2c32.pdf")
savefig("barbara_f2c32.png")

# Fig. 13d: eGHWT
eghwt32 = top_vectors_synthesis_2d(p, dvec_eghwt, loc_eghwt, GProws, GPcols, dmatrix)
savefig("barbara_eghwt32.pdf")
savefig("barbara_eghwt32.png")

#
# Zoom up the face region of those approximations (Fig. 14a)
#
row_zoom = 1:180
col_zoom = 330:450

# Haar
display(heatmap(haar32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_face_haar32.pdf")
savefig("barbara_face_haar32.png")

# Walsh
display(heatmap(walsh32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_face_walsh32.pdf")
savefig("barbara_face_walsh32.png")

# GHWT c2f
display(heatmap(c2f32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_face_c2f32.pdf")
savefig("barbara_face_c2f32.png")

# GHWT f2c
display(heatmap(f2c32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_face_f2c32.pdf")
savefig("barbara_face_f2c32.png")

# eGHWT
display(heatmap(eghwt32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_face_eghwt32.pdf")
savefig("barbara_face_eghwt32.png")

#
# Zoom up the left leg region of those approximations (Fig. 14b)
#
row_zoom = 300:512
col_zoom = 400:512

# Haar
display(heatmap(haar32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_lleg_haar32.pdf")
savefig("barbara_lleg_haar32.png")

# Walsh
display(heatmap(walsh32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_lleg_walsh32.pdf")
savefig("barbara_lleg_walsh32.png")

# GHWT c2f
display(heatmap(c2f32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_lleg_c2f32.pdf")
savefig("barbara_lleg_c2f32.png")

# GHWT f2c
display(heatmap(f2c32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_lleg_f2c32.pdf")
savefig("barbara_lleg_f2c32.png")

# eGHWT
display(heatmap(eghwt32[row_zoom, col_zoom],ratio=1, yaxis =:flip, showaxis = false, 
    ticks = false, colorbar = false, clims = (minimum(matrix), maximum(matrix)), color = :grays))
savefig("barbara_lleg_eghwt32.pdf")
savefig("barbara_lleg_eghwt32.png")

################################################################################
###################### Generate approximaation errors ##########################
################################################################################
# Generate Fig. 12b up to 50% of the coefficients retained.
DVEC =[ dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:] ]
T = [ "Graph Haar","Graph Walsh","GHWT_c2f", "GHWT_f2c", "eGHWT" ]
L = [ (:dashdot, :orange), (:dashdot, :blue), (:solid, :red), (:solid, :green),
      (:solid, :black)]
approx_error2(DVEC, T, L, 0.5)
savefig("barbara_approx_error.pdf")
savefig("barbara_approx_error.png")
