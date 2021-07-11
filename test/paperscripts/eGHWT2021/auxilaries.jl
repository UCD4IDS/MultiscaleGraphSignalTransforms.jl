################################################################################
####################### Approximation error plot ###############################
################################################################################
function approx_error2(DVEC::Array{Array{Float64,1},1}, 
                       T::Vector{String}, L::Vector{Tuple{Symbol,Symbol}},
                       frac::Float64 = 0.30)
    # This version plots the relative L2 errors against
    # the FRACTION of coefficients retained.
    plot(xaxis = "Fraction of Coefficients Retained", 
         yaxis = "Relative Approximation Error")
    for i = 1:length(DVEC)
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm
        er[er .== 0.0] .= minimum(er[er .!= 0.0]) # avoid blowup by taking log in the plot below
        # er is the relative L^2 error of the whole thing: length(er)=N.
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac), 
              label = T[i], line = L[i], linewidth = 2, grid = false)
    end
    display(current())
end


function approx_error3(DVEC::Array{Array{Float64,1},1},
                       T::Vector{String}, L::Vector{Tuple{Symbol,Symbol}},
                       N::Int64)
    # This version plots the relative L2 errors against
    # the NUMBER of coefficients retained.
    plot(xaxis = "Number of Coefficients Retained", 
         yaxis = "Relative Approximation Error")
    for i = 1:length(DVEC)
        dvec = DVEC[i]
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm 
        er[er .== 0.0] .= minimum(er[er .!= 0.0]) # avoid blowup by taking log in the plot below
        # er is the relative L^2 error of the whole thing: length(er)=N.
        plot!(0:N-1, er[1:N], yaxis=:log, label = T[i], line = L[i], 
              linewidth = 2, grid = false)
    end
    display(current())
end

################################################################################
############ Computing a weight matrix using the Gaussian affinity #############
################################################################################
function image_Gaussian_affinity(img::Matrix{Float64}, r::Int64, σ::Float64)
# Get the neighbors for affinity matrix computation
    l = 2*r + 1
    temp_x, temp_y = fill(1,l^2), fill(1,l^2)
    temp_xy = CartesianIndices((1:l,1:l))
    for i = 1:l^2
        temp_x[i], temp_y[i] = temp_xy[i][1], temp_xy[i][2]
    end
    # (r+1, r+1) is the index of the center location.
    temp_ind = ((temp_x .- (r + 1)).^2 + (temp_y .- (r + 1)).^2) .<= r^2
    # Now, temp_ind indicates those points withinin the circle of radius r 
    # from the center while neighbor_x, neighbor_y represent relative positions
    # of those points from the center.
    neighbor_x = temp_x[temp_ind] .- (r + 1) 
    neighbor_y = temp_y[temp_ind] .- (r + 1)
   # So, for any index (x, y), points within (x ± neighbor_x, y ± neighbor_y) are 
   # its neighbors for the purpose of calculating affinity

# Create affinity matrix
    m, n = size(img)
    sig = img[:]
    W = fill(0., (m*n,m*n))
    for i = 1:m*n
        cur = CartesianIndices((m,n))[i]
        for j = 1:length(neighbor_x) # assuming dim(neighbor_x) == dim(neighbor_y) 
            if 1 <= cur[1] + neighbor_x[j] <= m && 1 <= cur[2] + neighbor_y[j] <= n
                tempd = LinearIndices((m,n))[cur[1] + neighbor_x[j], cur[2] + neighbor_y[j]]
                W[i,tempd] = exp(-(sig[i] - sig[tempd])^2/σ)
            end
        end
    end
    return sparse(W)
end # end of the image_Gaussian_affinity function

################################################################################
######## Display top 9 basis vectors of various bases for an image data ########
################################################################################
function top_vectors_plot2(dvec::Array{Float64, 2}, m::Int64, n::Int64, BS::BasisSpec, GP::GraphPart;
                           clims::Tuple{Float64,Float64} = (-0.01, 0.01))

    # Get the indices from max to min in terms of absolute value of the coefs
    # Note that m*n == length(dvec); m and n are the image size.
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);

    # Set up the layout as 3 x 3 subplots
    plot(9, layout = (3,3), framestyle = :none, legend = false)

    # Do the display!
    for i=1:9
        dvecT = fill(0., size(dvec))
        dvecT[sorted_ind[i]] = 1
        f = ghwt_synthesis(dvecT, GP, BS)
        heatmap!(reshape(f, (m,n)), subplot=i, ratio=1, yaxis=:flip, 
                 showaxis=false, ticks = false, color = :grays, clims = clims)
    end
    display(current())
end
################################################################################

function top_vectors_plot3(dvec::Array{Float64, 2}, m::Int64, n::Int64, BS::BasisSpec,
    GP::GraphPart; clims::Tuple{Float64,Float64} = (-0.01, 0.01), K::Int64 = 9)

    # Get the indices from max to min in terms of absolute value of the coefs
    # Note that m*n == length(dvec); m and n are the image size.
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);

    # Set up the layout as 3 x 3 subplots
    plot(K, layout = (Int(sqrt(K)), Int(sqrt(K))), framestyle = :none, legend = false)

    # Do the display!
    for i=1:K
        dvecT = fill(0., size(dvec))
        dvecT[sorted_ind[i]] = 1
        f = ghwt_synthesis(dvecT, GP, BS)
        heatmap!(reshape(f, (m,n)), subplot=i, ratio=1, yaxis=:flip, 
                 showaxis=false, ticks = false, color = :grays, clims = clims)
    end
    display(current())
end
################################################################################

