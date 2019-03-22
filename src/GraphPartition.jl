module GraphPartition

using ..GraphSignal, SparseArrays, LinearAlgebra

include("partition_fiedler.jl")
include("utils.jl")

export GraphPart, partition_tree_fiedler, partition_tree_matrixDhillon

# GraphPart data structure and constructors
"""
    GP = GraphPart(ind::Vector{Int}, rs::Matrix{Int}, tag::Matrix{Int}, compinfo::Matrix{Int}, rsf2c::Matrix{Int}, tagf2c::Matrix{Int}, compinfof2c::Matrix{Int}, method::Vector{Symbol})

is a data structure for a GraphPart object containing the following fields:
* `ind::Vector{Int}`: ordering of the indices on the finest level
* `rs::Matrix{Int}`: regionstarInt (coarse-to-fine) <==> the index in `ind` of the first point in region number `i` is `rs[i]`
* `tag::Matrix{Int}`: tag info for the GHWT coarse-to-fine basis
* `compinfo::Matrix{Int}`: indicates whether the coefficient was formed from 2 coefficenInt (value is nonzero) or from only 1 coefficient (value is zero); when a scaling and Haar-like coefficient are formed, their corresponding values in compinfo indicate the number of nodes in each of the 2 subregions
* `rsf2c::Matrix{Int}`: the fine-to-coarse version of `rs`
* `tagf2c::Matrix{Int}`: the fine-to-coarse version of `tag`
* `compinfof2c::Matrix{Int}`: the fine-to-coarse version of `compinfo`
* `method::Symbol`: how the partition tree was constructed

The unsigned integer depends on the size of the underlying graph.

Copyright 2015 The RegenInt of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and modified by Naoki Saito, Feb. 7, 2017
Revised for two parameters by Naoki Saito, Feb. 24, 2017
"""
mutable struct GraphPart
    ind::Vector{Int}    # ordering of the indices on the finest level
    rs::Matrix{Int}     # `rs[i,j]` = the index in `ind` of the first
                       # point in Region `i` at level j
    tag::Matrix{Int}    # Int <: Int because max(tag) <= log2(max(T))
    compinfo::Matrix{Int}      # coef was from 2 coefs or from 1 coef
    rsf2c::Matrix{Int}         # f2c version of `rs`
    tagf2c::Matrix{Int}        # f2c version of `tag`
    compinfof2c::Matrix{Int}   # f2c version of `compinfo`
    method::Symbol           # specification of graph partition method

    # An inner constructor here.
    function GraphPart(ind::Vector{Int}, rs::Matrix{Int};
                       tag::Matrix{Int} = Matrix{Int}(undef, 0, 0),
                       compinfo::Matrix{Int} = Matrix{Int}(undef, 0, 0),
                       rsf2c::Matrix{Int} = Matrix{Int}(undef, 0, 0),
                       tagf2c::Matrix{Int} = Matrix{Int}(undef, 0, 0),
                       compinfof2c::Matrix{Int} = Matrix{Int}(undef, 0, 0),
                       method::Symbol = :unspecified)

        # Sanity checks
        if Base.size(rs, 1) != Base.length(ind) + 1
            @warn("size(rs,1) must be length(ind) + 1")
            @warn("both rs and ind now become null arrays!")
            ind = Vector{Int}(undef, 0)
            rs = Matrix{Int}(undef, 0, 0)
        end
        if Base.length(tag) != 0 && (Base.size(rs, 1) != Base.size(tag, 1) + 1
            || Base.size(rs, 2) != Base.size(tag, 2))
            @warn("tag size is inconsistent with rs size.")
            @warn("tag now becomes a null array!")
            tag = Matrix{Int}(undef, 0, 0)
        end
        if Base.length(compinfo) != 0 &&
            ( Base.size(rs, 1) != Base.size(compinfo, 1) + 1
             || Base.size(rs, 2) != Base.size(compinfo, 2) )
            @warn("compinfo size is inconsistent with rs size.")
            @warn("compinfo now becomes a null array!")
            compinfo = Matrix{Int}(undef, 0, 0)
        end
        if Base.length(rsf2c) != 0 && Base.length(rsf2c) != Base.length(rs)
            @warn("length(rsf2c) must be length(rs)")
            @warn("rsf2c now becomes a null array!")
            rsf2c = Matrix{Int}(undef, 0, 0)
        end
        if Base.length(tagf2c) != 0 && Base.length(tagf2c) != Base.length(tag)
            @warn("length(tagf2c) must be length(tag)")
            @warn("tagf2c now becomes a null array!")
            tagf2c = Matrix{Int}(undef, 0, 0)
        end
        if Base.length(compinfof2c) != 0 &&
            Base.length(compinfof2c) != Base.length(compinfo)
            @warn("length(compinfof2c) must be length(compinfo)")
            @warn("compf2c now becomes a null array!")
            compf2c = Matrix{Int}(undef, 0, 0)
        end
        new(ind, rs, tag, compinfo, rsf2c, tagf2c, compinfof2c, method)
    end # of an inner constructor GraphPart

 end # of type GraphPart


"""
    GP = partition_tree_fiedler(G::GraphSig, method::Symbol = :Lrw)

 Generate a partition tree for a graph using the Fiedler vector of either
 L (the unnormalized Laplacian) or L_rw (the random-walk normalized
 Laplacian).

### Input ArgumenInt
* `G::GraphSig`: an input GraphSig object
* `method::Symbol`: how the partition tree was constructed (default: :Lrw)

### Output Argument
* `GP::GraphPart`: an ouput GraphPart object

Copyright 2015 The RegenInt of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and modified by Naoki Saito, Feb. 7, 2017
"""
function partition_tree_fiedler(G::GraphSignal.GraphSig, method::Symbol = :Lrw)
    #
    # 0. Preliminary stuff
    #
    # constanInt
    N = G.length
    jmax = max(3 * floor(Int, log2(N)), 4) # jmax >= 4 is guaranteed.
    # This jmax is the upper bound of the true jmax; the true jmax
    # is computed as the number of columns of the matrix `rs` in the end.

    # TRACKING VARIABLES

    # `ind` records the way in which the nodes are indexed on each level
    ind = Vector{Int}(1:N)

    # `rs` stands for regionstarInt, meaning that the index in `ind` of the first
    # point in region number `i` is `rs[i]`
    rs = zeros(Int, N + 1, jmax)
    rs[1, :] .= 1
    rs[2, 1] = N + 1

    #
    # 1. Partition the graph to yield rs, ind, and jmax
    #
    j = 1                       # j=1 in julia <=> j=0 in theory
    regioncount = 0
    rs0 = 0                     # define here for the whole loops,
                                # which differs from MAIntAB.
    while regioncount < N
        #regioncount = countnz(rs[:,j]) - 1
        regioncount = count(!iszero,rs[:, j]) - 1 # the number of regions on level j
        if j == jmax  # add a column to rs for level j+1, if necessary
            rs = hcat(rs, vcat(Int(1), zeros(Int,N)))
            jmax = jmax + 1
        end
        # for tracking the child regions
        rr = 1
        for r = 1:regioncount   # cycle through the parent regions
            rs1 = rs[r, j]      # the start of the parent region
            rs2 = rs[r + 1, j] # 1 node after the end of the parent region
            n = rs2 - rs1   # the number of nodes in the parent region

            if n > 1            # regions with 2 or more nodes
                indrs = ind[rs1:(rs2 - 1)]
                # partition the current region
                (pm, ) = partition_fiedler(G.W[indrs,indrs], method = method)
                # determine the number of poinInt in child region 1
                n1 = sum(pm .> 0)
                # switch regions 1 and 2, if necessary, based on the sum of
                # edge weighInt to the previous child region (i.e., cousin)
                if r > 1
                    # MAIntAB: sum(sum(...)) instead of sum(...)
                    if sum(G.W[ind[rs0:(rs1 - 1)], indrs[pm .> 0]]) <
                        sum(G.W[ind[rs0:(rs1 - 1)], indrs[pm .< 0]])
                        pm = -pm
                        n1 = n - n1
                    end
                end
                # update the indexing
                ind[rs1:(rs1 + n1 - 1)] = indrs[pm .> 0]
                ind[(rs1 + n1):(rs2 - 1)] = indrs[pm .< 0]
                # update the region tracking
                rs[rr + 1, j + 1] = rs1 + n1
                rs[rr + 2, j + 1] = rs2
                rr = rr + 2
                rs0 = rs1 + n1
            elseif n == 1       # regions with 1 node
                rs[rr + 1, j + 1] = rs2
                rr = rr + 1
                rs0 = rs1
            end # of if n > 1 ... elseif n==1 construct
        end # of for r=1:regioncount
        j = j + 1
    end # of while regioncount < N statement

    #
    # 2. Postprocessing
    #
    # get rid of excess columns in rs
    rs = rs[:, 1:(j - 1)]         # in MAIntAB, it was rs(:,j:end) = [];
    # create a GraphPart object
    return GraphPart(ind, rs, method = method)
end # of function partition_tree_fiedler


"""
    GProws, GPcols = partition_tree_matrixDhillon(matrix::Matrix{Float64})

    Recursively partition the rows and columns of the matrix using the
    bipartite model proposed by Dhillon in "Co-clustering documents and words
    using Bipartite Spectral Graph Partitioning"

### Input Arguments
* `matrix::Matrix{Float64}`: an input matrix

### Output Argument
* `GProws::GraphPart`: partitioning using rows as samples
* `GPcols::GraphPart`: partitioning using cols as samples
"""
function partition_tree_matrixDhillon(matrix::Matrix{Float64})
  ## 0. Preliminary stuff

  # constants
  (rows, cols) = size(matrix)
  N = max(rows,cols)
  jmax = Int(max(3*floor(log2(N)),4))


  ### Tracking variables for the rows

  # order index
  ind_rows = Vector{Int}(1:rows)
  ind_cols = Vector{Int}(1:cols)

  # regionstarts
  rs_rows = zeros(Int, (N+1,jmax))
  rs_rows[1,:] .= 1
  rs_rows[2,1] = rows + 1

  rs_cols = zeros(Int, (N+1,jmax))
  rs_cols[1,:] .= 1
  rs_cols[2,1] = cols + 1


  ## 1. Partition the graph to yield rs, ind, and jmax
  j = 1
  regioncount_rows = 0
  regioncount_cols = 0
  while (regioncount_rows) < rows || (regioncount_cols < cols)
    # the number of true regions (i.e, has 1 or more nodes) on level j
    regioncount_rows = length(unique(vcat(rs_rows[:,j],[0]))) - 2
    regioncount_cols = length(unique(vcat(rs_cols[:,j],[0]))) - 2

    # the number of total regions (including those with no nodes)
    false_regioncount = max(count(!iszero, rs_rows[:,j]),count(!iszero, rs_cols[:,j])) - 1

    # add a column for level j+1, if necessary
    if j == jmax
      rs_rows = hcat(rs_rows,zeros(N+1))
      rs_cols = hcat(rs_cols,zeros(N+1))
      rs_rows[1,j+1] = 1
      rs_cols[1,j+1] = 1
      jmax = jmax + 1
    end

    # for tracking the child regions (for both rows and columns)
    rr  = 1

    # cycle through the parent regions
    for r = 1 : false_regioncount
      # the start of the parent region and 1 node after the end of the parent region
      rs1_rows = rs_rows[r,j]
      rs2_rows = rs_rows[r+1,j]
      rs1_cols = rs_cols[r,j]
      rs2_cols = rs_cols[r+1,j]

      # the number of nodes in the parent region
      n_rows = rs2_rows - rs1_rows
      n_cols = rs2_cols - rs1_cols

      if (n_rows <=1) && (n_cols > 1)
        indrs_rows = ind_rows[rs1_rows]
        indrs_cols = ind_cols[rs1_cols:(rs2_cols-1)]

        # compute the left and right singular vectors

        #(u_rows, v_cols) = second_largest_singular_vectors(reshape(matrix[indrs_rows,indrs_cols],(n_rows,n_cols)))


        # partition the current region
        pm_cols = partition_vector(matrix[indrs_rows,indrs_cols])

        # determine the number of points in child region 1
        n1_cols = sum(pm_cols .>0)

        # update the column indexing
        ind_cols[rs1_cols:(rs1_cols+n1_cols-1)] = indrs_cols[pm_cols.>0]
        ind_cols[(rs1_cols+n1_cols):(rs2_cols-1)] = indrs_cols[pm_cols.<0]

        # update the row region tracking
        rs_rows[rr+1,j+1] = rs1_rows
        rs_rows[rr+2,j+1] = rs2_rows

        # update the column region tracking
        rs_cols[rr+1, j+1] = rs1_cols + n1_cols
        rs_cols[rr+2, j+1] = rs2_cols

        rr = rr + 2

      elseif (n_rows > 1) && (n_cols <= 1)
        indrs_rows = ind_rows[rs1_rows:(rs2_rows-1)]
        indrs_cols = ind_cols[rs1_cols]

        # compute the left and right singular vectors
        #(u_rows, u_cols) = second_largest_singular_vectors(reshape(matrix[indrs_rows, indrs_cols],(n_rows,n_cols)))

        # partition the current region
        pm_rows = partition_vector(matrix[indrs_rows, indrs_cols])

        # determine the number of points in child region 1
        n1_rows = sum(pm_rows.>0)

        # update the row indexing
        ind_rows[rs1_rows:(rs1_rows+n1_rows-1)] = indrs_rows[pm_rows.>0]
        ind_rows[(rs1_rows+n1_rows):(rs2_rows-1)] = indrs_rows[pm_rows.<0]

        # update the row region tracking
        rs_rows[rr+1,j+1] = rs1_rows + n1_rows
        rs_rows[rr+2,j+1] = rs2_rows

        # update the column region tracking
        rs_cols[rr+1,j+1] = rs1_cols
        rs_cols[rr+2,j+1] = rs2_cols

        rr = rr + 2

      elseif (n_rows > 1) && (n_cols > 1)
        indrs_rows = ind_rows[rs1_rows:(rs2_rows-1)]
        indrs_cols = ind_cols[rs1_cols:(rs2_cols-1)]

        # compute the left and right singular vectors
        (u_rows,v_cols) = second_largest_singular_vectors(matrix[indrs_rows,indrs_cols])


        # partition the current region
        pm, = partition_fiedler_pm(vcat(u_rows, v_cols))
        pm_rows = pm[1:n_rows]
        pm_cols = pm[(n_rows+1):end]

        # determine the number of points in child region 1
        n1_rows = sum(pm_rows .> 0)
        n1_cols = sum(pm_cols .> 0)

        # update the row indexing
        ind_rows[rs1_rows:(rs1_rows + n1_rows - 1)] = indrs_rows[pm_rows .> 0]
        ind_rows[(rs1_rows + n1_rows):(rs2_rows - 1)] = indrs_rows[pm_rows .< 0]

        # update the column indexing
        ind_cols[rs1_cols:(rs1_cols + n1_cols - 1)] = indrs_cols[pm_cols .> 0]
        ind_cols[(rs1_cols + n1_cols):(rs2_cols - 1)] = indrs_cols[pm_cols .< 0]

        # update the row region tracking
        rs_rows[rr+1, j+1] = rs1_rows + n1_rows
        rs_rows[rr+2, j+1] = rs2_rows

        # update the column region tracking
        rs_cols[rr+1, j+1] = rs1_cols + n1_cols
        rs_cols[rr+2, j+1] = rs2_cols

        rr = rr + 2

      else
        rs_rows[rr+1, j+1] = rs2_rows
        rs_cols[rr+1, j+1] = rs2_cols
        rr = rr + 1
      end
      if rr > size(rs_rows,1) - 2
        rs_rows = vcat(rs_rows,zeros(Int,(2,size(rs_rows,2))))
      end
      if rr > size(rs_cols,1) - 2
        rs_cols = vcat(rs_cols,zeros(Int,(2,size(rs_cols,2))))
      end
    end

    j = j + 1
  end

  # get rid of excess columns in rs_rows and rs_cols
  rs_rows = rs_rows[:, 1:(j-1)]
  rs_cols = rs_cols[:, 1:(j-1)]
  jmax = size(rs_rows,2)

  #remove duplicates
  M_rows = size(rs_rows,1)
  M_cols = size(rs_cols,1)
  for j = jmax:-1:1
    ### ROWS
    # remove duplicate regionstarts
    unique_rows = unique(rs_rows[:,j])
    if unique_rows[1] == 0
      shift!(unique_rows)
    end
    rs_rows[:,j] = vcat(unique_rows, zeros(M_rows-length(unique_rows)))

    #remove duplicate columns of rs
    if (j < jmax) && (sum(rs_rows[:,j] .!= rs_rows[:,j+1]) == 0)
      rs_rows = rs_rows[:,setdiff(1:end,j+1)]
    end

    ### COLUMNS
    # remove duplicate regionscounts
    unique_cols = unique(rs_cols[:,j])
    if unique_cols[1] == 0
      shift!(unique_cols)
    end
    rs_cols[:,j] = vcat(unique_cols, zeros(M_cols-length(unique_cols)))

    # remove duplicate columns of rs
    if (j < jmax) && (sum(rs_cols[:,j].!=rs_cols[:,j+1]) == 0)
      rs_cols = rs_cols[:,setdiff(1:end,j+1)]
    end
  end

  # remove unnecessary rows and columns in the regionstarts
  if rows < M_rows - 1
    rs_rows = rs_rows[1:(rows+1),:]
  end
  if cols < M_cols - 1
    rs_cols = rs_cols[1:(cols+1),:]
  end

  # create GraphPart Objects
  GProws = GraphPart(ind_rows, rs_rows)
  GPcols = GraphPart(ind_cols, rs_cols)

  return GProws, GPcols

  # make sure the GraphPart objects meet our requirements
  #blablabla

end


function second_largest_singular_vectors(A::Matrix{Float64})
  # Compute the second largest left and right singular vectors of
  # An = D1^(-0.5)*A*D2^(-0.5)

  (rows, cols) = size(A)

  # compute D1 and D2 and make sure there are no zeros
  D1 = sum(A,dims = 2)
  D1[D1.==0] .= max(0.01, minimum(D1[D1.>0]/10))
  D2 = sum(A,dims = 1)
  D2[D2.==0] .= max(0.01, minimum(D2[D2.>0]/10))

  # compute the singular vectors
  (u,D,v) = svd(diagm(0 => D1[:].^(-0.5))*A*diagm(0 => D2[:].^(-0.5)), full = false)

  if (rows > 1) && (cols > 1)

  # extract the second left singular vector
  u = u[:,2]
  v = v[:,2]

# extract the 2nd singular vectors and multiply by D_i^(-0.5)
  u = D1[:].^(-0.5).*u
  v = D2[:].^(-0.5).*v

  return u, v


  end

end


function partition_vector(v:: Vector{Float64})
  # Partition vector into half, diveded by the mean of the values.
  l = size(v[:],1)
  m = mean(v)
  pm = zeros(l)
  pm[v.>m] .= 1
  pm[v.<=m] .= -1
  if size(unique(pm),1) == 1
    f = UInt(floor(l/2))
    pm[1:f] .= 1
    pm[f+1:end] .= -1
  end
  return pm
end


end # of module GraphPartition
