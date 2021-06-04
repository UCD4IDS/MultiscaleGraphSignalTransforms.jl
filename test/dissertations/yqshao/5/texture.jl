#Figure 5.9, 5.10
using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms, PyCall
py"""
import numpy
xs = numpy.load('../data/texture_mask.npy')
"""
xs = PyArray(py"xs"o)

matrix_gt = reshape(xs,(512,512))[1:4:512,1:4:512]
heatmap(matrix_gt,ratio=1, yaxis =:flip, axis = false, color = :grays, colorbar = false)
max1 = maximum(matrix_gt[:])
min1 = minimum(matrix_gt[:])
matrix_gt = (matrix_gt[:] .- min1)./(max1 - min1)

JLD2.@load "../data/handcut_images.jld2"
heatmap(tmp["block5"],ratio=1, yaxis =:flip, axis = false, color = :grays, colorbar = false)
matrix = tmp["block5"]

########################

sigma = 0.001
m, n = size(matrix)


############### Get the neighbors for affinity matrix computing
r = 5; # specify radius of neighbors
l = 2*r + 1

temp_x, temp_y = fill(1,l^2), fill(1,l^2)
temp_xy = CartesianIndices((1:l,1:l))
for i = 1:l^2
    temp_x[i], temp_y[i] = temp_xy[i][1], temp_xy[i][2]
end
temp_ind = ((temp_x .- (r + 1)).^2 + (temp_y .- (r + 1)).^2).<= r^2
neighbor_x = temp_x[temp_ind] .- (r + 1)
neighbor_y = temp_y[temp_ind] .- (r + 1)
# for any index(x,y), (x + neighbor_x, y - neigbor_y) are the neighbors to calculate affinity

################ Create affinity matrix
W = fill(0., (m*n,m*n))
for i = 1:m*n
    cur = CartesianIndices((m,n))[i]
    for j = 1:length(neighbor_x)
        if 1 <= cur[1] + neighbor_x[j] <= m && 1 <= cur[2] + neighbor_y[j] <= n
            tempd = LinearIndices((m,n))[cur[1] + neighbor_x[j], cur[2] + neighbor_y[j]]
            W[i,tempd] = exp(-(matrix_gt[i] - matrix_gt[tempd])^2/sigma)
        end
    end
end
W = sparse(W)

############### Preprocess to get G and GP
G = GraphSig(W, f = reshape(matrix,(length(matrix_gt),1)))
GP = partition_tree_fiedler(G)
dmatrix = ghwt_analysis!(G, GP=GP)

#####################
############# Haar
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

############# Walsh
BS_walsh = bs_walsh(GP)
dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)

############# GHWT_c2f
dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)

############# GHWT_f2c
dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)

############# eGHWT
dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)




function top_vectors_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart; clims::Tuple{Float64,Float64} = (-0.01,0.01))
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);
    plot(9,layout=(3,3),framestyle=:none, legend=false)
    for i=1:9
        dvecT = fill(0., size(dvec))
        #dvecT[sorted_ind[i]] = dvec_ghwt[sorted_ind[i]]
        dvecT[sorted_ind[i]] = 1
        f = ghwt_synthesis(dvecT, GP, BS)
        #print((maximum(f), minimum(f)))
        heatmap!(reshape(f, size(matrix)), subplot=i, ratio=1, yaxis=:flip, axis=false, color = :grays, clims = clims)
    end
    current()
end

#Figure 5.10b
top_vectors_plot(dvec_eghwt, BS_eghwt, GP)



function approx_error2(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0.3
    T = ["eGHWT Haar basis","eGHWT Walsh basis","GHWT_c2f", "GHWT_f2c", "eGHWT best basis"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :red),(:solid, :green),(:solid, :black)]
    for i = 1:5
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm # this is the relative L^2 error of the whole thing, i.e., its length is N
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac), label = T[i], line = L[i])
    end
end

#Figure 5.10a
approx_error2([dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:]])
