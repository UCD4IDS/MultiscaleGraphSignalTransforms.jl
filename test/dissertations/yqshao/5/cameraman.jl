#Figure 5.7, 5.8
using TestImages, Plots, SparseArrays, LinearAlgebra, Wavelets, MTSG
img = testimage("camera");
matrix = convert(Array{Float64,2}, img)[1:4:512,1:4:512]


heatmap(matrix,ratio=1, yaxis =:flip, axis = false, color = :grays, clim = (0,1), colorbar = false)

matrix_gt = matrix[:]

sigma = 0.007
#sigma = 0.07
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

############# Construct or search the specific basis
############# Haar
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

############# eGHWT
dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)


####################
function top_vectors_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart; clims::Tuple{Float64,Float64} = (-0.02,0.02))
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

#Figure 5.8a (re-run the code with sigma = 0.07 to generate figure 5.8b)
top_vectors_plot(dvec_eghwt, BS_eghwt, GP)

dvec_haar0007 = dvec_haar
dvec_eghwt0007 = dvec_eghwt
#re-run the code with sigma = 0.07 to generate dvec_haar007 and dvec_eghwt007
# dvec_haar007 = dvec_haar
# dvec_eghwt007 = dvec_eghwt


function approx_error2(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0.3
    T = ["Classical Haar transform", "eGHWT Haar basis (sigma = 0.07)", "eGHWT best basis (sigma = 0.07)", "eGHWT Haar basis (sigma = 0.007)", "eGHWT best basis (sigma = 0.007)"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :red),(:dashdot,:purple),(:solid,:black)]
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


dvec_classichaar = dwt(matrix, wavelet(WT.haar))

#Figure 5.7
approx_error2([dvec_classichaar[:], dvec_haar007[:], dvec_eghwt007[:], dvec_haar0007[:], dvec_eghwt0007[:]])
