using JLD2, SparseArrays, Images, ImageView, Gtk.ShortNames
@load "handcut_images.jld2"

#matrix = tmp["block5"]
#matrix_gt = tmp["block5_gt"][:]
#sigma = 0.0001

#matrix = tmp["block10"]
#matrix_gt = tmp["block10_gt"][:]
#sigma = 0.00001

matrix = tmp["lena"]
matrix_gt = tmp["lena_gt"][:]
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
            #if matrix_gt[i] == matrix_gt[tempd]
            #    W[i,tempd] = 1.0
            #end
            W[i,tempd] = exp(-(matrix_gt[i] - matrix_gt[tempd])^2/sigma)
        end
    end
end
W = sparse(W)
############### Preprocess to get G and GP
G = GraphSig(W, f = reshape(matrix,(length(matrix_gt),1)))
GP = partition_tree_fiedler(G)
dmatrix = ghwt_analysis!(G, GP=GP)

############### Visualization of cutting
BS_haar=bs_haar(GP)
dvec_haar=dmatrix2dvec(dmatrix, GP, BS_haar)

grid, frames, canvases = canvasgrid((4,4))

for i = 2:17
    dvecT = fill(0., size(dvec_haar))
    dvecT[i] = dvec_haar[i]
    f = ghwt_synthesis(dvecT, GP, BS_haar)
    f = reshape(f, size(matrix))
    imshow(canvases[i - 1], f)
end
win = Window(grid)
Gtk.showall(win)

############## Visuallization of Top ghwt vectors
dvec_ghwt,BS_ghwt = ghwt_bestbasis(dmatrix, GP,cfspec=1)
dvec_ghwt = reshape(dvec_ghwt, (length(dvec_ghwt),1))

sorted_ind = sortperm(abs.(dvec_ghwt[:]), rev = true);

grid, frames, canvases = canvasgrid((4,4))
for i = 1:16
    dvecT = fill(0., size(dvec_ghwt))
    dvecT[sorted_ind[i]] = dvec_ghwt[sorted_ind[i]]
    f = ghwt_synthesis(dvecT, GP, BS_ghwt)
    f = reshape(f, size(matrix))
    imshow(canvases[i], f)
end
win = Window(grid)
Gtk.showall(win)

############## Visualization of top ghwt vectors from tf selection

bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix[:,:,1], GP)
tf_ind = sortperm(abs.(bestbasis[:]), rev = true);

grid, frames, canvases = canvasgrid((4,4))
for i =1:16
    bestbasis_new = fill(0., size(bestbasis))
    bestbasis_new[tf_ind[i]] = bestbasis[tf_ind[i]]
    f, _ = tf_synthesis(bestbasis_new, bestbasis_tag, GP, G)
    f = reshape(f,(m,n))
    imshow(canvases[i], f)
end
win = Window(grid)
Gtk.showall(win)
