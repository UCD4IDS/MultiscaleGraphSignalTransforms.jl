"""
    pairclustering(𝚽::Matrix{Float64}, GP_star::GraphPart)

construct the GraphPart object of the primal graph via pair-clustering algorithm.

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `GP_star::GraphPart`: GraphPart object of the dual graph

# Output Argument
- `GP::GraphPart`: GraphPart object of the primal graph via pair-clustering

"""
function pairclustering(𝚽::Matrix{Float64}, GP_star::GraphPart)
    rs_dual = GP_star.rs
    inds_dual = GP_star.inds
    (N, jmax) = Base.size(inds_dual)

    # TRACKING VARIABLES

    # `ind` records the way in which the nodes are indexed on each level
    ind = Vector{Int}(1:N)

    # `inds` records the way in which the nodes are indexed on all levels
    inds = zeros(Int, N, jmax)
    inds[:, 1] = ind

    # `rs` stands for regionstart, meaning that the index in `ind` of the first
    # point in region number `i` is `rs[i]`
    rs = deepcopy(rs_dual)

    j = 1
    regioncount = 0
    while regioncount < N
        regioncount = count(!iszero, rs[:, j]) - 1
        if j == jmax  # add a column to rs for level j+1, if necessary
            rs = hcat(rs, vcat(Int(1), zeros(Int, N)))
            inds = hcat(inds, zeros(Int, N))
            jmax = jmax + 1
        end
        # for tracking the child regions
        rr = 1
        for r = 1:regioncount
            rs1 = rs[r, j]
            rs2 = rs[r + 1, j]
            n = rs2 - rs1
            if n > 1            # regions with 2 or more nodes
                indrs = ind[rs1:(rs2 - 1)]
                n1 = rs[rr + 1, j + 1] - rs1
                indrs_dual1 = inds_dual[rs1:(rs1 + n1 - 1), j + 1]
                indrs_dual2 = inds_dual[(rs1 + n1):(rs2 - 1), j + 1]
                # partition the current region
                (indp1, indp2) = partition_pc(𝚽, indrs, indrs_dual1, indrs_dual2)
                # update the indexing
                ind[rs1:(rs1 + n1 - 1)] = indrs[indp1]
                ind[(rs1 + n1):(rs2 - 1)] = indrs[indp2]
                # update the child region tracking
                rr += 2
            elseif n == 1
                rr += 1
            end
        end
        j = j + 1
        inds[:, j] = ind
    end

    # get rid of excess columns in rs
    rs = rs[:, 1:(j - 1)]
    # get rid of excess columns in inds
    inds = inds[:, 1:(j - 1)]

    return GraphPart(ind, rs; inds = inds)
end


"""
    partition_pc(𝚽::Matrix{Float64}, indrs::Vector{Int}, indrs_dual1::Vector{Int}, indrs_dual2::Vector{Int})

based on the partition on dual graph, partition the primal graph region
(indexed by `indrs`) into two pieces via the pair-clustering algorithm.

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `indrs::Vector{Int}`: indices of the primal graph region to be partitioned
- `indrs_dual1::Vector{Int}`: index1 of the dual graph region's partition result
- `indrs_dual2::Vector{Int}`: index2 of the dual graph region's partition result

# Output Argument
- `indp1::Vector{Int}`: index1 of `indrs` for the partition result
- `indp2::Vector{Int}`: index2 of `indrs` for the partition result

"""
function partition_pc(𝚽::Matrix{Float64}, indrs::Vector{Int}, indrs_dual1::Vector{Int}, indrs_dual2::Vector{Int})
    n = length(indrs)
    n1 = length(indrs_dual1)
    E = 𝚽[indrs, :].^2
    score = sum(E[:, indrs_dual1], dims = 2) - sum(E[:, indrs_dual2], dims = 2)
    indp1 = sortperm(score[:]; rev = true)[1:n1]
    indp2 = setdiff(1:n, indp1)
    return indp1, indp2
end


"""
    pc_ngwp(𝚽::Matrix{Float64}, GP_star::GraphPart, GP::GraphPart)

construct pair-clustering NGWP and GP.tag in place.

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `GP_star::GraphPart`: GraphPart object of the dual graph
- `GP::GraphPart`: GraphPart object of the primal graph

# Output Argument
- `wavelet_packet::Array{Float64,3}`: the pair-clustering NGWP. The first
    index is for selecting wavelets at a fixed level; the second index is for
    selecting the level `j`; the third index is for selecting elements in the
    wavelet vector.

"""
function pc_ngwp(𝚽::Matrix{Float64}, GP_star::GraphPart, GP::GraphPart)
    rs = GP_star.rs
    inds_dual = GP_star.inds
    inds_primal = GP.inds
    (N, jmax) = Base.size(inds_dual)

    GP_star.tag = zeros(Int, N, jmax)
    GP_star.tag[:, 1] = Vector{Int}(0:(N - 1))

    wavelet_packet = zeros(N, jmax, N)
    wavelet_packet[:, 1, :] = Matrix{Float64}(I, N, N)

    for j = 2:jmax
        regioncount = count(!iszero, rs[:, j]) - 1
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            GP_star.tag[indr, j] = Vector{Int}(0:(length(indr) - 1))
            wavelet_packet[indr, j, :] = const_proj_wavelets(𝚽, inds_primal[indr, j], inds_dual[indr, j])'
        end
    end

    return wavelet_packet
end
