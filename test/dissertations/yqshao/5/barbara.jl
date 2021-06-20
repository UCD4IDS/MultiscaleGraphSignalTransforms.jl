#Figure 5.4, 5.5
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms

JLD2.@load "../data/spie_data.jld2"
matrix = vars["barbara"]

#face
row_zoom = 1:180
col_zoom = 330:450

#right leg
row_zoom = 300:512
col_zoom = 400:512

#full image
row_zoom = 1:512
col_zoom = 1:512

heatmap(matrix[row_zoom, col_zoom],ratio=1, yaxis =:flip, axis = false, color = :grays)
#savefig("original.pdf")


### initialize the regular balanced binary partition and compute the expanding coefficients
dmatrix, GProws, GPcols = ghwt_tf_init_2d_Lindberg(matrix)

### Compute the coefficients of different basis we will compare
############# Haar
BS_haar_rows = bs_haar(GProws)
BS_haar_cols = bs_haar(GPcols)
dvec_haar, loc_haar = BS2loc(dmatrix, GProws, GPcols, BS_haar_rows, BS_haar_cols)

############# Walsh
BS_walsh_rows = bs_walsh(GProws)
BS_walsh_cols = bs_walsh(GPcols)
dvec_walsh, loc_walsh = BS2loc(dmatrix, GProws, GPcols, BS_walsh_rows, BS_walsh_cols)

############# GHWT (i.e., regular haar-walsh wavelet packet dictionary)
fcols, jmax_col = size(GProws.tag);
frows, jmax_row = size(GPcols.tag);
dmatrix_rows = reshape(dmatrix, (frows, jmax_row, fcols*jmax_col))
dmatrix_cols = Array{Float64,3}(reshape(dmatrix',(fcols, jmax_col, frows*jmax_row)))
############# c2f bestbasis
dvec_c2f_rows, BS_c2f_rows = ghwt_c2f_bestbasis(dmatrix_rows, GProws)
dvec_c2f_cols, BS_c2f_cols = ghwt_c2f_bestbasis(dmatrix_cols, GPcols)
dvec_c2f, loc_c2f = BS2loc(dmatrix, GProws, GPcols, BS_c2f_rows, BS_c2f_cols)

############# f2c bestbasis
dvec_f2c_rows, BS_f2c_rows = ghwt_f2c_bestbasis(dmatrix_rows, GProws)
dvec_f2c_cols, BS_f2c_cols = ghwt_f2c_bestbasis(dmatrix_cols, GPcols)
dvec_f2c, loc_f2c = BS2loc(dmatrix, GProws, GPcols, BS_f2c_rows, BS_f2c_cols)

#############

##################### tf bestbasis
dvec_tf, loc_tf = ghwt_tf_bestbasis_2d(matrix, GProws, GPcols)

################################################################################
##################### Visualize the results of synthesis########################
################################################################################


### function to plot the image synthesized by top p vectors
function top_vectors_synthesis_2d(p::Int64, dvec::Vector{Float64}, loc::Matrix{Int64}, GProws::GraphPart, GPcols::GraphPart, dmatrix::Matrix{Float64})
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0

    matrix_syn = ghwt_tf_synthesis_2d(dvecT, loc, GProws, GPcols)
    heatmap(matrix_syn[row_zoom, col_zoom],ratio=1, yaxis =:flip, axis = false, color = :grays)
    mse = norm(matrix - matrix_syn,2)^2/length(matrix)
    psnr = -10*log10(mse)
    return psnr
end

################################################################################
percent = 1/32;
p = Int64(floor(percent*length(matrix)))

### haar
# Figure 5.4(a)
top_vectors_synthesis_2d(p, dvec_haar, loc_haar, GProws, GPcols, dmatrix)
#savefig("synthesis_haar_1_32.pdf")

### walsh
top_vectors_synthesis_2d(p, dvec_walsh, loc_walsh, GProws, GPcols, dmatrix)
#savefig("synthesis_walsh_1_32.pdf")

### c2f
# Figure 5.4(b)
top_vectors_synthesis_2d(p, dvec_c2f, loc_c2f, GProws, GPcols, dmatrix)
#savefig("synthesis_c2f_1_32.pdf")

### f2c
# Figure 5.4(c)
top_vectors_synthesis_2d(p, dvec_f2c, loc_f2c, GProws, GPcols, dmatrix)
#savefig("synthesis_f2c_1_32.pdf")

### tf
# Figure 5.4(d)
top_vectors_synthesis_2d(p, dvec_tf, loc_tf, GProws, GPcols, dmatrix)
#savefig("synthesis_tf_1_32.pdf")

# To generate Figure 5.5, just change the row_zoom, and col_zoom at the head of this file.

################################################################################
####################### Approximation error plot################################
################################################################################
function approx_error(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0:0.01:0.3
    T = ["Haar","Walsh","GHWT_c2f", "GHWT_f2c", "eGHWT"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :red),(:solid, :green),(:solid, :black)]
    for i = 1:5
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2, rev = true)

        er = fill(0., length(frac))

        for j = 1:length(frac)
            p = Int64(floor(frac[j]*N))
            er[j] = sqrt(dvec_norm^2 - sum(dvec_sort[1:p]))/dvec_norm
        end

        plot!(frac, er, yaxis=:log, xlims = (0.,0.3), label = T[i], line = L[i])
    end
end

approx_error([dvec_haar, dvec_walsh, dvec_c2f, dvec_f2c, dvec_tf])
current()



#################
### Generate the relative l2 error when approximating by 1/32 coefficients
################
DVEC = [dvec_haar, dvec_walsh, dvec_c2f, dvec_f2c, dvec_tf];
for i = 1:5
    dvec = DVEC[i]
    N = length(dvec)
    dvec_norm = norm(dvec,2)
    dvec_sort = sort(dvec.^2, rev = true)
    p = Int64(floor(N./32))
    print(sqrt(dvec_norm^2 - sum(dvec_sort[1:p]))/dvec_norm)
end
