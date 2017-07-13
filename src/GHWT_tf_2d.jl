include("PartitionTreeMatrixDhillon.jl")
include("partition_fiedler.jl")

function tf_core_2d(dmatrix::Matrix{Float64}, tag::Matrix{<:Any}, tag_r::Matrix{<:Any})
  # the core iteration of the time-frequency analysis algorithm
  costfun = function(x) return sum(abs.(x),1) end
  (m,n) = size(tag)
  (p,q) = size(dmatrix)
  dmatrix_new = zeros(p,q)
  tag_new = -1*ones(m,n-1)
  tag_r_new = -1*ones(m,n-1)
  tag_tf = zeros(p,q) # indicating if time or frequency domain is chosen

  # the i-th element of "totallengths" indicate the number of total useful entries
  # from level 1 to level i-1 in tag matrix
  totallengths = zeros(Int,n)
  totallengths[1] = 0
  for i=2:n
    totallengths[i] = sum(tag[:,i-1].!=-1) + totallengths[i-1]
  end


  s = 1
  for j=1:(n-1)
    region = unique(tag_r[:,j])
    region = region[region.!=-1]
    regioncount = length(region)
    i = 1
    for r = 1:regioncount
      parent_ind = find(tag_r[:,j].==region[r])
      child1_ind = find(tag_r[:,j+1].==region[r]*2)
      child2_ind = find(tag_r[:,j+1].==(region[r]*2+1))
      for k = parent_ind[1]:parent_ind[end]
        if tag[k,j]%2 == 0
          tempind1 = find(tag[child1_ind,j+1].==tag[k,j]/2) #index of child1
          tempind2 = find(tag[child2_ind,j+1].==tag[k,j]/2) #index of child2
          if isempty(tempind2)
            timecos = costfun(dmatrix[child1_ind[1]-1+tempind1+totallengths[j+1],:])
          elseif isempty(tempind1)
            timecos = costfun(dmatrix[child2_ind[1]-1+tempind2+totallengths[j+1],:])
          else
            timecos = costfun(vcat(dmatrix[child1_ind[1]-1+tempind1+totallengths[j+1],:], dmatrix[child2_ind[1]-1+tempind2+totallengths[j+1],:]))
          end # children correspond to time direction

          if k == parent_ind[end]
            freqcos = costfun(dmatrix[k+totallengths[j],:]')
          elseif tag[k+1,j] == tag[k,j] + 1 #check if it is paired with some odd tag
            freqcos = costfun(vcat(dmatrix[k+totallengths[j],:]', dmatrix[k+1+totallengths[j],:]'))
          else
            freqcos = costfun(dmatrix[k+totallengths[j],:]')
          end # parents correspond to frequency direction

          dmatrix_new[s,:], tag_tf[s,:] = findmin(vcat(freqcos, timecos),1)
          tag_tf[s,:] = tag_tf[s,:].%2
          tag_new[i,j] = tag[k,j]/2
          tag_r_new[i,j] = tag_r[k,j]
          i = i+1
          s = s+1
        end
      end
    end
  end

  #remove the redundant rows
  ind = (sum(tag_new.!=-1,2).!=0)[:]
  tag_new = tag_new[ind,:]
  tag_r_new = tag_r_new[ind,:]

  dmatrix_new = dmatrix_new[1:s-1,:]
  tag_tf = tag_tf[1:s-1,:]

  return dmatrix_new, tag_new, tag_r_new, tag_tf

end



"""
    dmatrix, GProws, GPcols = ghwt_tf_init_2d(matrix::Matrix{Float64})

    Partition matrix first to get GProws and GPcols. Then expand matrix
    in two directions to get dmatrix

### Input Arguments
* `matrix::Matrix{Float64}`: an input matrix

### Output Argument
* `GProws::GraphPart`: partitioning using rows as samples
* `GPcols::GraphPart`: partitioning using cols as samples
* `dmatrix::matrix{Float64}`: expansion coefficients of matrix
"""
function ghwt_tf_init_2d(matrix::Matrix{Float64})
  #Partition matrix recuisively using Dhillon's method
  GProws, GPcols = PartitionTreeMatrixDhillon(matrix)
  ghwt_core!(GProws)
  ghwt_core!(GPcols)


  matrix_re = matrix[GProws.ind, GPcols.ind]

  (N, jmax_row) = size(GProws.rs)
  N = N-1

  # expand on the row direction
  (frows, fcols) = size(matrix_re)
  dmatrix = zeros((N,jmax_row,fcols))
  dmatrix[:,jmax_row,:] = matrix_re

  ghwt_core!(GProws, dmatrix)

  dmatrix = reshape(dmatrix, (frows*jmax_row, fcols))'

  # expand on the column direction
  (N, jmax_col) = size(GPcols.rs)
  N = N - 1
  dmatrix2 = zeros(N,jmax_col,size(dmatrix,2))
  dmatrix2[:, jmax_col, :] = dmatrix
  ghwt_core!(GPcols, dmatrix2)

  dmatrix = reshape(dmatrix2, (fcols*jmax_col, frows*jmax_row))'

  return dmatrix, GProws, GPcols
end



"""
    Bbasis, infovec = ghwt_tf_bestbasis_2d(dmatrix::Matrix{Float64}, GProws::GraphPart, GPcols::GraphPart)

    Find the best-basis with time-frequency analysis.

### Input Arguments
* `dmatrix::Matrix{Float64}`: the expansion coefficients matrix
* `GProws::GraphPart`: partitioning using rows as samples
* `GPcols::GraphPart`: partitioning using cols as samples

### Output Argument
* `Bbasis::matrix{Float64}`: the value of bestbasis
* `infovec::matrix{Float64}`: the location of bestbasis in dmatrix

"""
function ghwt_tf_bestbasis_2d(dmatrix::Matrix{Float64}, GProws::GraphPart, GPcols::GraphPart)

  tag_r_row = rs_to_region(GProws.rs, GProws.tag)
  tag_r_col = rs_to_region(GPcols.rs, GPcols.tag)

  jmax_col = size(GPcols.tag,2)
  jmax_row = size(GProws.tag,2)

  # tag information using columns as samples
  TAG_col = Vector(jmax_col)
  TAG_col[1] = GPcols.tag
  TAG_r_col = Vector(jmax_col)
  TAG_r_col[1] = tag_r_col

  # tag information using rows as samples
  TAG_row = Vector(jmax_row)
  TAG_row[1] = GProws.tag
  TAG_r_row = Vector(jmax_row)
  TAG_r_row[1] = tag_r_row

  # recording if time or frequency domain is chosen
  TAG_tf = Matrix(jmax_row,jmax_col)

  # storing the cost of the basis
  DMATRIX = Matrix(jmax_row, jmax_col)

  DMATRIX[1,1] = dmatrix

  # comparison in the row direction
  for i=1:(jmax_row - 1)
    DMATRIX[i+1,1], TAG_row[i+1], TAG_r_row[i+1], TAG_tf[i+1,1] = tf_core_2d(DMATRIX[i,1], TAG_row[i], TAG_r_row[i])
    TAG_tf[i+1,1] = Matrix{Int8}(TAG_tf[i+1,1] + 2)
  end

  # comparison in the column direction
  for i=1:(jmax_col - 1)
    DMATRIX[1,i+1], TAG_col[i+1], TAG_r_col[i+1], TAG_tf[1,i+1] = tf_core_2d(DMATRIX[1,i]',TAG_col[i], TAG_r_col[i])
    DMATRIX[1,i+1] = DMATRIX[1,i+1]'
    TAG_tf[1,i+1] = Matrix{Int8}(TAG_tf[1,i+1]')
  end

  DMATRIX[1,1] = 0 # release the memory

  # comparison in both directions
  for i = 2:jmax_row
    for j = 2:jmax_col
      (temp_col,_,_,temp_tf_col) = tf_core_2d(DMATRIX[i,j-1]', TAG_col[j-1], TAG_r_col[j-1])
      temp_col = temp_col'
      temp_tf_col = temp_tf_col'
      (temp_row,_,_,temp_tf_row) = tf_core_2d(DMATRIX[i-1,j], TAG_row[i-1], TAG_r_row[i-1])

      (p,q) = size(temp_row)
      temp_tf_col = temp_tf_col[:]
      temp_tf_row = temp_tf_row[:]
      DMATRIX[i,j], TAG_tf[i,j] = findmin(vcat(temp_col[:]', temp_row[:]'),1)
      DMATRIX[i,j] = reshape(DMATRIX[i,j], (p,q))
      ind_r = ((TAG_tf[i,j].%2).==0)[:]
      ind_c = ((TAG_tf[i,j].%2).==1)[:]

      TAG_tf[i,j][ind_r] = temp_tf_row[ind_r]+2


      TAG_tf[i,j][ind_c] = temp_tf_col[ind_c]

      TAG_tf[i,j] = Matrix{Int8}(reshape(TAG_tf[i,j], (p,q)))

      DMATRIX[i-1,j] = 0 # release the memory

    end

    DMATRIX[i,1] = 0 # release the memory

  end

  for j = 2: jmax_col
    DMATRIX[jmax_row, j] = 0 # release the memory
  end


  #recover the bestbasis
  infovec = Array{Int}([jmax_row; jmax_col; 1; 1])

  for iter = 1:(jmax_row + jmax_col - 2)
    newinfovec1 = -1*ones(4,size(infovec,2))
    newinfovec2 = -1*ones(4,size(infovec,2))
    for h = 1:size(infovec,2)
      i = infovec[1,h]
      j = infovec[2,h]
      p = infovec[3,h]
      q = infovec[4,h]
      # if the bestbasis comes from column direction
      if TAG_tf[i,j][p,q] == 0 || TAG_tf[i,j][p,q] == 1
        s1 = 0
        s2 = 0
        l = 0
        k = 0

        for s = 1 : size(TAG_col[j],2)
          s2 = sum(TAG_col[j][:,s].!=-1) + s1
          if (s1 < q) && (q <= s2)
            l = s
            k = q - s1
            break
          else
            s1 = deepcopy(s2)
          end
        end
        region = TAG_r_col[j][k,l]
        tag = TAG_col[j][k,l]
        if TAG_tf[i,j][p,q] == 1 #frequency domain was chosen
          k1 = find((TAG_col[j-1][:,l].==tag*2).&(TAG_r_col[j-1][:,l].==region))
          k2 = find((TAG_col[j-1][:,l].==tag*2+1).&(TAG_r_col[j-1][:,l].==region))
        end
        if TAG_tf[i,j][p,q] == 0 #time domain was chosen
          k1 = find((TAG_col[j-1][:,l+1].==tag).&(TAG_r_col[j-1][:,l+1].==region*2))
          k2 = find((TAG_col[j-1][:,l+1].==tag).&(TAG_r_col[j-1][:,l+1].==region*2+1))
          l = l + 1
        end
        s1 = 0
        for s = 1:(l-1)
          s1 = s1 + sum(TAG_col[j-1][:,s].!=-1)
        end
        q1 = s1 + k1
        q2 = s1 + k2
        if isempty(q2)
          newinfovec1[:,h] = [i;j-1;p;q1]
        elseif isempty(q1)
          newinfovec2[:,h] = [i;j-1;p;q2]
        else
          newinfovec1[:,h] = [i;j-1;p;q1]
          newinfovec2[:,h] = [i;j-1;p;q2]
        end
      end

      # if the bestbasis comes from the row direction
      if TAG_tf[i,j][p,q] == 2 || TAG_tf[i,j][p,q] == 3
        s1 = 0
        s2 = 0
        l = 0
        k = 0


        for s = 1:size(TAG_row[i],2)
          s2 = sum(TAG_row[i][:,s].!=-1) + s1
          if s1 < p && p <= s2
            l = s
            k = p - s1
            break
          else
            s1 = deepcopy(s2)
          end
        end
        region = TAG_r_row[i][k,l]
        tag = TAG_row[i][k,l]
        if TAG_tf[i,j][p,q] == 3 #frequency domain was chosen
          k1 = find((TAG_row[i-1][:,l].==tag*2).&(TAG_r_row[i-1][:,l].==region))
          k2 = find((TAG_row[i-1][:,l].==tag*2+1).&(TAG_r_row[i-1][:,l].==region))
        end
        if TAG_tf[i,j][p,q] == 2 #time domain was chosen
          k1 = find((TAG_row[i-1][:,l+1].==tag).&(TAG_r_row[i-1][:,l+1].==region*2))
          k2 = find((TAG_row[i-1][:,l+1].==tag).&(TAG_r_row[i-1][:,l+1].==region*2+1))
          l = l + 1
        end
        s1 = 0
        for s = 1:(l-1)
          s1 = s1 + sum(TAG_row[i-1][:,s].!=-1)
        end
        p1 = s1 + k1
        p2 = s1 + k2
        if isempty(p2)
          newinfovec1[:,h] = [i-1;j;p1;q]
        elseif isempty(p1)
          newinfovec2[:,h] = [i-1;j;p2;q]
        else
          newinfovec1[:,h] = [i-1;j;p1;q]
          newinfovec2[:,h] = [i-1;j;p2;q]
        end
      end
    end
    infovec=Matrix{Int}(hcat(newinfovec1[:,newinfovec1[1,:].!=-1], newinfovec2[:,newinfovec2[1,:].!=-1]))

  end

  #obtain the values of the bestbasis coefficients
  Bbasis = dmatrix[sub2ind(size(dmatrix),infovec[3,:],infovec[4,:])]'

  return Bbasis, vcat(infovec[[3],:],infovec[[4],:])
end
