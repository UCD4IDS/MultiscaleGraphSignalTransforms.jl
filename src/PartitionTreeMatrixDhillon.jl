include("utils.jl")


"""
    GProws, GPcols = PartitionTreeMatrixDhillon(matrix::Matrix{Float64})

    Recursively partition the rows and columns of the matrix using the
    bipartite model proposed by Dhillon in "Co-clustering documents and words
    using Bipartite Spectral Graph Partitioning"

### Input Arguments
* `matrix::Matrix{Float64}`: an input matrix

### Output Argument
* `GProws::GraphPart`: partitioning using rows as samples
* `GPcols::GraphPart`: partitioning using cols as samples
"""
function PartitionTreeMatrixDhillon(matrix::Matrix{Float64})
  ## 0. Preliminary stuff

  # constants
  (rows, cols) = size(matrix)
  N = max(rows,cols)
  Tl = ind_class(N)
  jmax = Int(max(3*floor(log2(N)),4))


  ### Tracking variables for the rows

  # order index
  ind_rows = Vector{Tl}(1:rows)
  ind_cols = Vector{Tl}(1:cols)

  # regionstarts
  rs_rows = zeros(Tl, (N+1,jmax))
  rs_rows[1,:] = 1
  rs_rows[2,1] = rows + 1

  rs_cols = zeros(Tl, (N+1,jmax))
  rs_cols[1,:] = 1
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
    false_regioncount = max(countnz(rs_rows[:,j]),countnz(rs_cols[:,j])) - 1

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
        rs_rows = vcat(rs_rows,zeros(Tl,(2,size(rs_rows,2))))
      end
      if rr > size(rs_cols,1) - 2
        rs_cols = vcat(rs_cols,zeros(Tl,(2,size(rs_cols,2))))
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
  Ts = tag_class(max(size(rs_rows,2),size(rs_cols,2)))
  GProws = GraphPart{Tl,Ts}(ind_rows, rs_rows)
  GPcols = GraphPart{Tl,Ts}(ind_cols, rs_cols)

  return GProws, GPcols

  # make sure the GraphPart objects meet our requirements
  #blablabla

end


function second_largest_singular_vectors(A::Matrix{Float64})
  # Compute the second largest left and right singular vectors of
  # An = D1^(-0.5)*A*D2^(-0.5)

  (rows, cols) = size(A)

  # compute D1 and D2 and make sure there are no zeros
  D1 = sum(A,2)
  D1[D1.==0] = max(0.01, minimum(D1[D1.>0]/10))
  D2 = sum(A,1)
  D2[D2.==0] = max(0.01, minimum(D2[D2.>0]/10))

  # compute the singular vectors
  (u,D,v) = svd(diagm(D1[:].^(-0.5))*A*diagm(D2[:].^(-0.5)), thin = true)

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
  pm[v.>m] = 1
  pm[v.<=m] = -1
  if size(unique(pm),1) == 1
    f = UInt(floor(l/2))
    pm[1:f] = 1
    pm[f+1:end] = -1
  end
  return pm
end
