using Plots, SparseArrays, JLD2, LinearAlgebra, Wavelets, MTSG

JLD2.@load "../4.2.1/presentation.jld2"
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

#face_smaller
row_zoom = 71:170
col_zoom = 341:440

heatmap(matrix[row_zoom, col_zoom],ratio=1, yaxis =:flip, axis = false, color = :grays)
#savefig("original.pdf")
matrix = matrix[row_zoom, col_zoom]

### initialize the regular balanced binary partition and compute the expanding coefficients
GPcols = PartitionTreeMatrix_unbalanced_haar_binarytree(matrix, 9, 3.)
GProws = PartitionTreeMatrix_unbalanced_haar_binarytree(Matrix{Float64}(transpose(matrix)), 8, 3.)

dmatrix = ghwt_analysis_2d(matrix, GProws, GPcols)

### Compute the coefficients of different basis we will compare
############# Haar
BS_haar_rows = bs_haar(GProws)
BS_haar_cols = bs_haar(GPcols)
dvec_haar, loc_haar = BS2loc(dmatrix, GProws, GPcols, BS_haar_rows, BS_haar_cols)

# ############# Walsh
# BS_walsh_rows = bs_walsh(GProws)
# BS_walsh_cols = bs_walsh(GPcols)
# dvec_walsh, loc_walsh = BS2loc(dmatrix, GProws, GPcols, BS_walsh_rows, BS_walsh_cols)
#
# ############# GHWT (i.e., regular haar-walsh wavelet packet dictionary)
# fcols, jmax_col = size(GProws.tag);
# frows, jmax_row = size(GPcols.tag);
# dmatrix_rows = reshape(dmatrix, (frows, jmax_row, fcols*jmax_col))
# dmatrix_cols = Array{Float64,3}(reshape(dmatrix',(fcols, jmax_col, frows*jmax_row)))
# ############# c2f bestbasis
# dvec_c2f_rows, BS_c2f_rows = ghwt_c2f_bestbasis(dmatrix_rows, GProws)
# dvec_c2f_cols, BS_c2f_cols = ghwt_c2f_bestbasis(dmatrix_cols, GPcols)
# dvec_c2f, loc_c2f = BS2loc(dmatrix, GProws, GPcols, BS_c2f_rows, BS_c2f_cols)
#
# ############# f2c bestbasis
# dvec_f2c_rows, BS_f2c_rows = ghwt_f2c_bestbasis(dmatrix_rows, GProws)
# dvec_f2c_cols, BS_f2c_cols = ghwt_f2c_bestbasis(dmatrix_cols, GPcols)
# dvec_f2c, loc_f2c = BS2loc(dmatrix, GProws, GPcols, BS_f2c_rows, BS_f2c_cols)

#############

##################### tf bestbasis
dvec_tf, loc_tf = ghwt_tf_bestbasis_2d(matrix, GProws, GPcols)


##################### classical haar
dvec_classichaar = dwt(matrix, wavelet(WT.haar))

################################################################################
#########################
function approx_error(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0:0.01:0.3
    T = ["Classical Haar transform", "eGHWT Haar basis", "eGHWT best basis"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :black)]
    for i = 1:3
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
        print(er[26])
    end
end

################################################################################
####################### Approximation error plot################################
################################################################################
### function to plot the approximation error curve
function approx_error2(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0.3
    T = ["Classical Haar transform", "eGHWT Haar basis", "eGHWT best basis"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :black)]
    for i = 1:3
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm # this is the relative L^2 error of the whole thing, i.e., its length is N
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac), label = T[i], line = L[i])
    end
end

### Figure 7
#approx_error([dvec_classichaar[:], dvec_haar[:], dvec_tf[:]])
approx_error2([dvec_classichaar[:], dvec_haar[:], dvec_tf[:]])
current()
savefig("figure7.pdf")
