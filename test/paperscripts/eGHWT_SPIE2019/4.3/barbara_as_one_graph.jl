using Plots, SparseArrays, JLD2, LinearAlgebra, Wavelets, MultiscaleGraphSignalTransforms

#img = load("test\\8.jpg")
#heatmap(img,ratio=1, yaxis =:flip, axis = false, color = :gray)

JLD2.@load "../spie_data.jld2"
matrix = vars["barbara"][1:4:512,1:4:512]

heatmap(matrix,ratio=1, yaxis =:flip, axis = false, color = :grays, clim = (0,1))

matrix_gt = matrix[:]

sigma = 0.01
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

G_tmp = Glevel(G,GP,4)
heatmap(reshape(G_tmp.f[:],(128,128)))
############# Construct or search the specific basis
############# Haar
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

############# eGHWT
dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)

############# Walsh
BS_walsh = bs_walsh(GP)
dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)

############# GHWT_c2f
dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)

############# GHWT_f2c
dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)





###################################################################################
############## Visuallization of Top basis vector (multiple)#######################
###################################################################################
### visualization fo the first 12 vectors
### function to plot top 12 vectors
function top_vectors_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart; clims::Tuple{Float64,Float64} = (-0.02,0.02))
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);
    plot(12,layout=(3,4),framestyle=:none, legend=false)
    for i=1:12
        dvecT = fill(0., size(dvec))
        #dvecT[sorted_ind[i]] = dvec_ghwt[sorted_ind[i]]
        dvecT[sorted_ind[i]] = 1
        f = ghwt_synthesis(dvecT, GP, BS)
        #print((maximum(f), minimum(f)))
        heatmap!(reshape(f, size(matrix)), subplot=i, ratio=1, yaxis=:flip, axis=false, color = :balance, clims = clims)
    end
    current()
end

### haar
top_vectors_plot(dvec_haar, BS_haar, GP)
### walsh
top_vectors_plot(dvec_walsh, BS_walsh, GP)
### c2f
top_vectors_plot(dvec_c2f, BS_c2f, GP)
### f2c
top_vectors_plot(dvec_f2c, BS_f2c, GP)
### eghwt
top_vectors_plot(dvec_eghwt, BS_eghwt, GP)

###################################################################################
############## Visuallization of Top basis vector (single)#########################
###################################################################################
### function to plot the first i-th vector
### Note that here only one vector is plotted instead of multiple vectors.
function top_vector_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart, G::GraphSig, i::Int64; clims::Tuple{Float64,Float64} = (-0.02,0.02))
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);
    dvecT = fill(0., size(dvec))
    #dvecT[sorted_ind[i]] = dvec_ghwt[sorted_ind[i]]
    dvecT[sorted_ind[i]] = 1
    f, GS = ghwt_synthesis(dvecT, GP, BS, G)
    heatmap(reshape(f, size(matrix)), ratio=1, yaxis=:flip, axis=false, color = :grays, clims = clims)
end
########################## plot top i-th vector
i = 2
### haar
top_vector_plot(dvec_haar, BS_haar, GP, G, i)
### walsh
top_vector_plot(dvec_walsh, BS_walsh, GP, G, i)
### c2f
top_vector_plot(dvec_c2f, BS_c2f, GP, G, i)
### f2c
top_vector_plot(dvec_f2c, BS_f2c, GP, G, i)
### eghwt
top_vector_plot(dvec_eghwt, BS_eghwt, GP, G, i)





#######################################################################################
### generate figure used for spie paper
##########################################################################################
# Figure 8b
sorted_ind = sortperm(abs.(dvec_eghwt[:]), rev = true)
plot(4, layout=(2,2), framestyle=:none, legend=false)

for i = 1:4
    dvec_T = fill(0., size(dvec_eghwt))
    dvec_T[sorted_ind[i+1]] = 1
    f = ghwt_synthesis(dvec_T, GP, BS_eghwt)
    heatmap!(reshape(f, size(matrix)), subplot=i, ratio=1, yaxis=:flip, axis=false, color = :grays, clims = (-0.01,0.01))
end
current()

savefig("barbara_top4.pdf")


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

#########################
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

dvec_classichaar = dwt(matrix, wavelet(WT.haar))

### Figure 8(a)
#approx_error([dvec_classichaar[:], dvec_haar[:], dvec_eghwt[:]])
approx_error2([dvec_classichaar[:], dvec_haar[:], dvec_eghwt[:]])
current()
