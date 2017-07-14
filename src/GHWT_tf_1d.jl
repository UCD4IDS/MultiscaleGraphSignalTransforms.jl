module GHWT_tf_1d

include("utils.jl")

using ..GraphSignal, ..GraphPartition, ..BasisSpecification

include("common.jl")

export tf_core, ghwt_tf_bestbasis, tf_basisrecover, tf_threshold, tf_synthesis


"""
    dmatrix_new, tag_new, tag_r_new,tag_tf = tf_core(dmatrix::Matrix{Float64},tag::Matrix{<:Any},tag_r::Matrix{<:Any})

Perform one iteration of the time-frequancy analysis

### Input Arguments
* `dvec::Matrix{Float64}`: the expansion coefficients corresponding to the chosen basis or sums of cost of them
* `tag::Matrix{<:Any}`: indicating the tag of the coefficient of the element
* `tag_r::Matrix{<:Any}`: indicating the region of the elements

### Output Arguments
* `dmatrix_new::Matrix{Float64}`: Sums of cost, each of which is computed from at most four elements in dmatrix. The
smaller of the cost from the frequency direction or time direction is chosen.
* `tag_new::Matrix{<:Any}`: indicating the tag of the coefficient of the element
* `tag_r_new::Matrix{<:Any}`: indicating the region of the elements
* `tag_tf::Matrix{<:Any}`: indicating whether time cost (0) or frequency cost (1) is chosen
"""
function  tf_core(dmatrix::Matrix{Float64},tag::Matrix{<:Any},tag_r::Matrix{<:Any})
  costfun = cost_functional(1)
  (m,n) = size(dmatrix)
  dmatrix_new = zeros(m,n-1)
  tag_new = -1*ones(m,n-1)
  tag_r_new = -1*ones(m,n-1)
  tag_tf = -1*ones(m,n-1)

  for j = 1:(n-1)  # j is the level
    region = unique(tag_r[:,j]) # indicate the region number on level j
    deleteat!(region,find(region.==-1)) #remove complementary number -1
    regioncount = length(region)
    i = 1
    for r = 1:regioncount # iterate region wise
      parent_ind = find(tag_r[:,j] .== region[r])
      child1_ind = find(tag_r[:,j+1] .== region[r]*2)
      child2_ind = find(tag_r[:,j+1] .== region[r]*2+1)
      for k = parent_ind[1]:parent_ind[end]
        if tag[k,j]%2 == 0 # only look at even tags, since all the odds are paired with evens

          tempind1 = find(tag[child1_ind,j+1].==tag[k,j]/2) #index of child1
          tempind2 = find(tag[child2_ind,j+1].==tag[k,j]/2) #index of child2
          if isempty(tempind2)
            timecos = costfun(dmatrix[child1_ind[1]-1+tempind1,j+1])
          elseif isempty(tempind1)
            timecos = costfun(dmatrix[child2_ind[1]-1+tempind2,j+1])
          else
            timecos = costfun([dmatrix[child1_ind[1]-1+tempind1,j+1] dmatrix[child2_ind[1]-1+tempind2,j+1]])
          end # children correspond to time direction

          if k == parent_ind[end]
            freqcos = costfun(dmatrix[k,j])
          elseif tag[k+1,j] == tag[k,j] + 1 #check if it is paired with some odd tag
            freqcos = costfun([dmatrix[k,j] dmatrix[k+1,j]])
          else
            freqcos = costfun(dmatrix[k,j])
          end # parents correspond to frequency direction

          if timecos <= freqcos
            dmatrix_new[i,j] = timecos
            tag_tf[i,j] = 0
          else
            dmatrix_new[i,j] = freqcos
            tag_tf[i,j] = 1
          end
          tag_new[i,j] = tag[k,j]/2
          tag_r_new[i,j] = tag_r[k,j]
          i = i+1
        end
      end
    end
  end

  ind = sum(tag_new.!=-1,2).!=0 #the index of rows not completely equal to -1
  ind = ind[:,1]
  dmatrix_new = dmatrix_new[ind,:]
  tag_new = Matrix{Int32}(tag_new[ind,:])
  tag_r_new = Matrix{Int32}(tag_r_new[ind,:])
  tag_tf = Matrix{Int8}(tag_tf[ind,:])
  return dmatrix_new, tag_new, tag_r_new,tag_tf
end




"""
    bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix, GP)
Find the best basis with time-frequency analysis

### Input Arguments
* `dmatrix::Matrix{Float64}`: the expansion coefficients corresponding to the chosen basis
* `GP::GraphPart`: an input GraphPart object

### Output Arguments
* `bestbasis::Matrix{Float64}`: Sums of cost, each of which is computed from at most four elements in dmatrix. The
smaller of the cost from the frequency direction or time direction is chosen.
* `bestbasis_tag::Matrix{Int8}`: Indicating which of the coefficients are chosen
"""

function ghwt_tf_bestbasis(dmatrix, GP)
  tag = GP.tag
  tag_r = rs_to_region(GP.rs, GP.tag) # From rs and tag, find the region information.
  (m,n) = size(dmatrix)
  TAG = Vector(n) # store the tag information in every iteration
  TAG_r = Vector(n) # store the tag_r information in every iteration
  TAG_tf = Vector(n) # store the tag_tf information in every iteration
  DMATRIX = Vector(n) # store the cost in every iteration

  TAG[1] = tag
  TAG_r[1] = tag_r
  TAG_tf[1] = zeros(m,n)
  DMATRIX[1] = dmatrix[:,:,1]

  for i=1:(n-1) # compare the costs by iterations
    DMATRIX[i+1], TAG[i+1], TAG_r[i+1],TAG_tf[i+1] = tf_core(DMATRIX[i],TAG[i],TAG_r[i])
  end


  bestbasis_tag = TAG_tf[n]
  for i = (n-1):-1:1 # recover the bestbasis from tag_tf
    bestbasis_tag = tf_basisrecover(TAG_tf[i], bestbasis_tag, TAG[i], TAG_r[i])
  end

  bestbasis = deepcopy(dmatrix)
  bestbasis[bestbasis_tag.==-1] = 0
  bestbasis_tag = Matrix{Int8}( bestbasis_tag + 1 )
  return bestbasis, bestbasis_tag
end


"""
tag_tf_b_new = tf_basisrecover(tag_tf_b::Matrix{<:Any}, tag_tf_f::Matrix{<:Any}, tag::Matrix{<:Any}, tag_r::Matrix{<:Any})
One iteration in recovering the bestbasis in ghwt_tf_bestbasis method.

### Input Arguments
* `tag_new::Matrix{<:Any}`: indicating the tag of the coefficient of the element
* `tag_r_new::Matrix{<:Any}`: indicating the region of the elements
* `tag_tf_f::Matrix{<:Any}`:in the current iteration, indicating whether time cost (0) or frequency cost (1) is chosen
* `tag_tf_b::Matrix{<:Any}`:in the previous iteration, indicating whether time cost (0) or frequency cost (1) is chosen


### Output Arguments
* `tag_tf_b_new::Matrix{Float64}`: same as tag_tf_b, but only the costs corresponding to bestbasis is chosen
"""


function tf_basisrecover(tag_tf_b::Matrix{<:Any}, tag_tf_f::Matrix{<:Any}, tag::Matrix{<:Any}, tag_r::Matrix{<:Any})
  m,n = size(tag_tf_b)
  tag_tf_b_new = -1*ones(m,n)

  for j = 1:(n-1)
    region = unique(tag_r[:,j])
    deleteat!(region,find(region.==-1))
    regioncount = length(region)
    i = 1
    for r = 1:regioncount
      parent_ind = find(tag_r[:,j] .== region[r])
      child1_ind = find(tag_r[:,j+1] .== region[r]*2)
      child2_ind = find(tag_r[:,j+1] .== region[r]*2+1)
      for k = parent_ind[1]:parent_ind[end]
        if tag[k,j]%2 == 0
          if tag_tf_f[i,j] == 0
            tempind1 = find(tag[child1_ind,j+1].==tag[k,j]/2)
            tempind2 = find(tag[child2_ind,j+1].==tag[k,j]/2)
            if isempty(tempind2)
              tag_tf_b_new[child1_ind[1]-1+tempind1,j+1] = tag_tf_b[child1_ind[1]-1+tempind1,j+1]
            elseif isempty(tempind1)
              tag_tf_b_new[child2_ind[1]-1+tempind2,j+1] = tag_tf_b[child2_ind[1]-1+tempind2,j+1]
            else
              tag_tf_b_new[child1_ind[1]-1+tempind1,j+1] = tag_tf_b[child1_ind[1]-1+tempind1,j+1]
              tag_tf_b_new[child2_ind[1]-1+tempind2,j+1] = tag_tf_b[child2_ind[1]-1+tempind2,j+1]
            end
          elseif tag_tf_f[i,j] == 1
            if k == parent_ind[end]
              tag_tf_b_new[k,j] = tag_tf_b[k,j]
            elseif tag[k+1,j] == tag[k,j]+1
              tag_tf_b_new[k,j] = tag_tf_b[k,j]
              tag_tf_b_new[k+1,j] = tag_tf_b[k+1,j]
            else
              tag_tf_b_new[k,j] = tag_tf_b[k,j]
            end
          end
          i = i+1
        end
      end
    end
  end
  return tag_tf_b_new
end



"""
 bestbasis_new = tf_threshold(bestbasis::Matrix{Float64}, GP::GraphPart, keep::Float64, SORH::String)

 Thresholding the coefficients of bestbasis.

### Input Arguments
 *   `bestbasis::Matrix{Float64}`        the matrix of expansion coefficients
 *   `SORH::String`        use soft ('s') or hard ('h') thresholding
 *   `keep::Float64`        a fraction between 0 and 1 which says how many coefficients should be kept
 *   `GP::GraphPart`          a GraphPart object, used to identify scaling coefficients

### Output Argument
 *   `bestbasis_new::Matrix{Float64}`       the thresholded expansion coefficients
"""

function tf_threshold(bestbasis::Matrix{Float64}, GP::GraphPart, keep::Float64, SORH::String)

  tag = GP.tag
  if keep > 1 || keep < 0
    error("keep should be floating point between 0~1")
  end
  kept = UInt32(round(keep*size(bestbasis,1)))
  dvec_S = sort(abs.(bestbasis[:]), rev = true)
  T = dvec_S[kept + 1]
  bestbasis_new = deepcopy(bestbasis[:])
  indp = bestbasis_new.> T           #index for coefficients > T
  indn = bestbasis_new.< -1*T        #index for coefficients < -T

  # hard thresholding
  if SORH == "h" || SORH == "hard"
    bestbasis_new[.~(indp .| indn)] = 0

  # soft thresholding
  elseif SORH == "s" || SORH == "soft"
    bestbasis_new[(.~indp) .& (.~indn) .& (tag[:].!=0)] = 0
    bestbasis_new[indp .& (tag[:].!=0)] = bestbasis[indp .& (tag[:].!=0)] - T
    bestbasis_new[indn .& (tag[:].!=0)] = bestbasis[indn .& (tag[:].!=0)] + T
  end

  bestbasis_new = reshape(bestbasis_new,size(bestbasis))
end




"""
    (f, GS) = tf_synthesis(bestbasis::Matrix{Float64},bestbasis_tag::Matrix{<:Any},GP::GraphPart,G::GraphSig)

Given a vector of GHWT expansion coefficients and info about the graph
partitioning and the choice of basis, reconstruct the signal

### Input Arguments
* `bestbasis::Matrix{Float64}`: the expansion coefficients corresponding to the chosen basis
* 'bestbasis_tag::Matrix{<:Any}': the location of the best basis coefficients in bestbasis matrix
* `GP::GraphPart`: an input GraphPart object
* `G::GraphSig`: an input GraphSig object

### Output Arguments
* `f::Matrix{Float64}`: the reconstructed signal(s)
* `GS::GraphSig`: the reconstructed GraphSig object
"""
function tf_synthesis(bestbasis::Matrix{Float64},bestbasis_tag::Matrix{<:Any},GP::GraphPart,G::GraphSig)
  tag = GP.tag
  rs = GP.rs
  bestbasis_new = deepcopy(bestbasis)
  jmax = size(rs,2)
  for j = 1:(jmax-1)
    regioncount = countnz(rs[:,j]) - 1
    for r = 1:regioncount
      # the index that marks the start of the first subregion
      rs1 = rs[r,j]

      # the index that is one after the end of the second subregion
      rs3 = rs[r+1,j]

      # the number of points in the current region
      n = rs3 - rs1

      # only proceed forward if the coefficients do not exist
      if countnz(bestbasis_tag[rs1:(rs3-1),j]) !=0
        if n == 1
          # scaling coefficient
          if bestbasis_tag[rs1,j] == 1 # check ind
            bestbasis_new[rs1,j+1] = bestbasis_new[rs1,j]
            bestbasis_tag[rs1,j+1] = 1
          end
        elseif n > 1
          # the index that marks the start of the second subregion
          rs2 = rs1 + 1
          while rs2 < rs3 && tag[rs2, j+1] != 0
            rs2 = rs2 +1
          end

          # only one child
          if rs2 == rs3
            if bestbasis_tag[rs1:rs3-1,j] == 1
              bestbasis_new[rs1:rs3-1,j+1] = bestbasis_new[rs1:rs3-1,j]
              bestbasis_tag[rs1:rs3-1,j+1] = 1
            end

          else

            # the number of points in the first subregion
            n1 = rs2-rs1
            # the number of points in the second subregion
            n2 = rs3-rs2

            # scaling coefficients
            if bestbasis_tag[rs1,j] == 1 && bestbasis_tag[rs1+1,j] == 1 # check if it is the coefficients of best basis
              bestbasis_new[rs1,j+1] = (sqrt(n1) *bestbasis_new[rs1,j] + sqrt(n2)*bestbasis_new[rs1+1,j])/sqrt(n)
              bestbasis_new[rs2,j+1] = (sqrt(n2) *bestbasis_new[rs1,j] - sqrt(n1)*bestbasis_new[rs1+1,j])/sqrt(n)
              bestbasis_tag[rs1,j+1] = 1
              bestbasis_tag[rs2,j+1] = 1
            end

            ### HAAR-LIKE & WALSH-LIKE coefficients

            # search through the remaining coefficients in each subregion
            parent = rs1 + 2
            child1 = rs1 + 1
            child2 = rs2 + 1
            while child1 < rs2 || child2 < rs3
              # subregion 1 has the smaller tag
              if child2 == rs3 || (tag[child1,j+1] < tag[child2, j+1] && child1 < rs2)
                if bestbasis_tag[parent,j]==1 # check if it is the coefficients of best basis
                  bestbasis_new[child1,j+1] = bestbasis_new[parent,j]
                  bestbasis_tag[child1, j+1] =1
                end
                child1 = child1 + 1
                parent = parent + 1

              # subregion 2 has the smaller tag
              elseif child1 == rs2 || (tag[child2, j+1] < tag[child1, j+1] && child2 < rs3)
                if bestbasis_tag[parent, j] == 1 # check if it is the coefficients of best basis
                  bestbasis_new[child2, j+1] = bestbasis_new[parent, j]
                  bestbasis_tag[child2, j+1] = 1
                end
                child2 = child2 + 1
                parent = parent + 1

                # both subregions have the same tag
              else
                if bestbasis_tag[parent,j] == 1 && bestbasis_tag[parent+1, j] == 1 # check if it is the coefficients of best basis
                  bestbasis_new[child1,j+1] = (bestbasis_new[parent,j] + bestbasis_new[parent+1,j])/sqrt(2)
                  bestbasis_new[child2,j+1] = (bestbasis_new[parent,j] - bestbasis_new[parent+1,j])/sqrt(2)
                  bestbasis_tag[child1,j+1]=1
                  bestbasis_tag[child2,j+1]=1
                end
                child1 = child1 + 1
                child2 = child2 + 1
                parent = parent + 2
              end
            end
          end
        end
      end
    end
  end
  ftemp = bestbasis_new[:,end]

  f = zeros(size(ftemp))
  f[GP.ind] = ftemp
  f = reshape(f,(length(f),1)) # reorder f
  GS = deepcopy(G)
  replace_data!(GS, f) # create the new graph signal
  return f, GS
end


end
