module GHWT_tf_1d

include("utils.jl")

using ..GraphSignal, ..GraphPartition, ..BasisSpecification, LinearAlgebra

include("common.jl")

export ghwt_tf_bestbasis, tf_threshold, tf_synthesis


"""
    coeffdict = tf_init(dmatrix::Matrix{Float64},GP::GraphPart)

Store the expanding coeffcients from matrix into a list of dictionary (inbuilt hashmap in Julia)

### Input Arguments
* `dmatrix`: The expanding GHWT coefficients of all levels corresponding to input GP
* `GP::GraphPart`: an input GraphPart object

### Output Arguments
* `coeffdict`: The expanding GHWT coeffcients stored in a list of "dictionary" (inbuilt hashmap in Julia),
* `coeffdict`: The entry `coeffdict[j][(k,l)]` corresponds to the coefficient of basis-vector on level j with region k and tag l.

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function tf_init(dmatrix::Matrix{Float64},GP::GraphPart)
    # Obtain the tag information
    tag = convert(Array{Int64,2},GP.tag)

    # Obtain the region information
    tag_r = convert(Array{Int64,2},rs_to_region(GP.rs, GP.tag))

    (m,n) = size(dmatrix)

    # Initialize coeffdict
    coeffdict = Array{Dict{Tuple{Int64,Int64},Float64}}(undef, n)

    # Fill in the values with the rule that `coeffdict[j][(k,l)]` represents
    # the coefficient of basis-vector on level j with region k and tag l.
    for i = 1:n
        coeffdict[i]= Dict{Tuple{Int64,Int64},Float64}((tag_r[j,i],tag[j,i]) => dmatrix[j,i] for j=1:m)
    end
    return coeffdict
end





"""
    coeffdict_new,tag_tf = tf_core_new(coeffdict::Array{Dict{Tuple{Int64,Int64},Float64},1})


One forward iteration of time-frequency adapted GHWT method. For each entry in `coeffdict_new`, we compare
two (or one) entries in 'coeffdict' on time-direction and two (or one) entries in 'coeffdict' on frequency-direction.
Those two groups reprensent the same subspace. We compare the cost-functional value of them and choose the smaller one
as a new entry in 'coeffdict_new'.


### Input Arguments
* `coeffdict`: The entries of which reprensents the cost functional value of some basis-vectors' coefficients.

### Output Arguments
* `coeffdict_new`: The entries of which represents the cost functional value of some basis-vectors' coefficients
* `tag_tf`: Indicating whether the time-direction (0) or frequency direction (1) was chosen for each entry in coeffdict_new.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function tf_core_new(coeffdict::Array{Dict{Tuple{Int64,Int64},Float64},1})
    # We choose '1-norm' as the optional cost functional
    costfun = cost_functional(1)
    jmax = length(coeffdict)

    # Initialization
    coeffdict_new = Array{Dict{Tuple{Int64,Int64},Float64}}(undef, jmax-1)
    tag_tf = Array{Dict{Tuple{Int64,Int64},Bool}}(undef, jmax-1)

    # Iterate through levels
    for j = 1:(jmax-1)
        # the temporary dictionary, which will be the j-th level of `coeffdict_new`
        temp_coeff = Dict{Tuple{Int64,Int64},Float64}()

        # the temporary dictionary, which will be the j-th level of 'tag_tf'
        temp_tf = Dict{Tuple{Int64,Int64},Bool}()

        # iterate through the entries in coeffdict on level j
        for key in keys(coeffdict[j])
            # coeffdict[j][k,l] represent the entry on level j with region k and tag l
            k = key[1] # region index
            l = key[2] # tag index

            # only look at the entry with even tag l to avoid duplication
            if l%2 == 0

                # search for the (j,k,l+1) entry.
                # the (j,k,l) and (j,k,l+1) span the same subspace as (j+1,2*k,l/2) and (j+1,2*k+1,l/2)
                # (j,k,l) and (j,k,l+1) are `frequency-direction`
                # (j,k,l) and (j,k,l+1) are `time-direction`

                if haskey(coeffdict[j],(k,l+1)) # check for degenerate case ((j,k,l+1) doesn't exist)
                    freqcos = costfun([coeffdict[j][(k,l)],coeffdict[j][(k,l+1)]])
                else
                    freqcos = costfun([coeffdict[j][(k,l)]])
                end

                if ~haskey(coeffdict[j+1],(2*k,l/2)) # check for degenerate case ((j+1,2*k,l/2) or (j+1,2*k+1,l/2) doesn't exist)
                    timecos = costfun([coeffdict[j+1][(2*k+1,l/2)]])
                elseif ~haskey(coeffdict[j+1],(2*k+1,l/2))
                    timecos = costfun([coeffdict[j+1][(2*k,l/2)]])
                else
                    timecos = costfun([coeffdict[j+1][(2*k+1,l/2)],coeffdict[j+1][(2*k,l/2)]])
                end

                # compare the cost-functional value and record into 'tag_tf'
                if timecos <= freqcos
                    temp_coeff[(k,l/2)] = timecos
                    temp_tf[(k,l/2)] = false
                else
                    temp_coeff[(k,l/2)] = freqcos
                    temp_tf[(k,l/2)] = true
                end
            end
        end
        coeffdict_new[j] = temp_coeff
        tag_tf[j] = temp_tf
    end
    return coeffdict_new,tag_tf
end






"""
    tag_tf_b_new = tf_basisrecover_new(tag_tf_b::Array{Dict{Tuple{Int64,Int64},Bool}},tag_tf_f::Array{Dict{Tuple{Int64,Int64},Bool}})


One backward iteration of time-frequency adapted GHWT method to recover the best-basis from the `tag_tf`s recorded.

### Input Arguments
* `tag_tf_b`: The `dictionary` recording the time-or-frequency information on some iteration 'i' in the main algorithm
* `tag_tf_f`: The `dictionary` recording the time-or-frequency information on some iteration 'i+1' in the main algorithm


### Output Arguments
* `tag_tf_b_new`: The updated 'tag_tf_b'. Eventually the 'tag_tf' on iteration 1 will represent the selected best-basis


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function tf_basisrecover_new(tag_tf_b::Array{Dict{Tuple{Int64,Int64},Bool}},tag_tf_f::Array{Dict{Tuple{Int64,Int64},Bool}})
    # Initialization
    jmax = length(tag_tf_b)
    #tag_tf_b_new = Array{Dict{Tuple{Int64,Int64},Bool}}(jmax)
    tag_tf_b_new = Array{Dict{Tuple{Int64,Int64},Bool}}(undef, jmax)
    for j = 1:jmax
        tag_tf_b_new[j] = Dict{Tuple{Int64,Int64},Bool}()
    end

    # Iterate on the levels
    for j = 1:(jmax-1)
        for key in keys(tag_tf_f[j])
            k = key[1]
            l = key[2]

            # The entries on frequency-direction are selected
            if tag_tf_f[j][(k,l)] == true
                if ~haskey(tag_tf_b[j],(k,2*l))
                    tag_tf_b_new[j][(k,2*l+1)] =  tag_tf_b[j][(k,2*l+1)]
                elseif  ~haskey(tag_tf_b[j],(k,2*l+1))
                    tag_tf_b_new[j][(k,2*l)] = tag_tf_b[j][(k,2*l)]
                else
                    tag_tf_b_new[j][(k,2*l)] =  tag_tf_b[j][(k,2*l)]
                    tag_tf_b_new[j][(k,2*l+1)] =  tag_tf_b[j][(k,2*l+1)]
                end
            else
            # The entries on time-direction are selected
                if ~haskey(tag_tf_b[j+1],(2*k,l))
                    tag_tf_b_new[j+1][(2*k+1,l)] =  tag_tf_b[j+1][(2*k+1,l)]
                elseif  ~haskey(tag_tf_b[j+1],(2*k+1,l))
                    tag_tf_b_new[j+1][(2*k,l)] =  tag_tf_b[j+1][(2*k,l)]
                else
                    tag_tf_b_new[j+1][(2*k+1,l)] =  tag_tf_b[j+1][(2*k+1,l)]
                    tag_tf_b_new[j+1][(2*k,l)] =  tag_tf_b[j+1][(2*k,l)]
                end
            end
        end
    end
    return tag_tf_b_new
end


"""
bestbasis_tag_matrix, bestbasis = ghwt_tf_bestbasis_new(dmatrix::Matrix{Float64},GP::GraphPart)

Implementation of time-frequency adapted GHWT method.
Modified from the algorithm in paper 'A Fast Algorithm for Adapted Time Frequency Tilings' by Christoph M Thiele and Lars F Villemoes.

### Input Arguments
### Input Arguments
* `dmatrix`: The expanding GHWT coefficients of all levels corresponding to input GP.
* `GP::GraphPart`: an input GraphPart object.

### Output Arguments
* `bestbasis_tag_matrix`: binary 0-1 matrix indicating the location of best-basis in dmatrix
* `bestbasis`: same size as dmatrix, but only coefficients of best-basis vectors are nonzero

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function ghwt_tf_bestbasis(dmatrix::Matrix{Float64},GP::GraphPart)
    (m,n) = size(dmatrix)

    # Initialization. Store the expanding coeffcients from matrix into a list of dictionary (inbuilt hashmap in Julia)
    # The entry `coeffdict[j][(k,l)]` corresponds to the coefficient of basis-vector on level j with region k and tag l.
    tag = convert(Array{Int64,2},GP.tag)
    tag_r = convert(Array{Int64,2},rs_to_region(GP.rs, GP.tag))
    (m,n) = size(dmatrix)
    coeffdict = Array{Dict{Tuple{Int64,Int64},Float64}}(undef, n)
    for i = 1:n
        coeffdict[i]= Dict{Tuple{Int64,Int64},Float64}((tag_r[j,i],tag[j,i]) => dmatrix[j,i] for j=1:m)
    end

    # TAG_tf stores the time-or-frequency information on every iteration
    # COEFFDICT stores the corresponding cost-functional values.
    COEFFDICT = Vector(undef, n)
    TAG_tf = Vector(undef, n)

    # Initialization of the first iteration
    COEFFDICT[1] = coeffdict
    TAG_tf_init = Array{Dict{Tuple{Int64,Int64},Bool}}(undef, n)
    for j = 1:n
        TAG_tf_init[j] = Dict(key => true for key in keys(coeffdict[j]))
    end
    TAG_tf[1] = TAG_tf_init

    # Iterate forward. For each entry in `COEFFDICT[i+1]`, we compare two (or one) entries in 'COEFFDICT[i]' on time-direction and two (or one) on frequency-direction.
    # Those two groups reprensent the same subspace. We compare the cost-functional value of them and choose the smaller one as a new entry in 'COEFFDICT[i+1]'
    # 'TAG_tf[i+1]' records if frequency-direction (1) or time-direction (0) was chosen.
    for i = 1:(n-1)
        COEFFDICT[i+1], TAG_tf[i+1] = tf_core_new(COEFFDICT[i])
    end

    # Iterate backward with the existing tag_tf information to recover the best-basis.
    bestbasis_tag = TAG_tf[n]
    for i = (n-1):-1:1
        bestbasis_tag = tf_basisrecover_new(TAG_tf[i],bestbasis_tag)
    end

    # Change the data structure from dictionary to matrix
    bestbasis = zeros(m,n)
    bestbasis_tag_matrix = zeros(m,n)
    for j = 1:n
        for i = 1:m
            k = tag_r[i,j]
            l = tag[i,j]
            if haskey(bestbasis_tag[j],(k,l))
                bestbasis_tag_matrix[i,j] = 1
                bestbasis[i,j] = dmatrix[i,j]
            end
        end
    end
    return bestbasis, bestbasis_tag_matrix
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
    regioncount = count(!iszero, rs[:,j]) - 1
    for r = 1:regioncount
      # the index that marks the start of the first subregion
      rs1 = rs[r,j]

      # the index that is one after the end of the second subregion
      rs3 = rs[r+1,j]

      # the number of points in the current region
      n = rs3 - rs1

      # only proceed forward if the coefficients do not exist
      if count(!iszero, bestbasis_tag[rs1:(rs3-1),j]) !=0
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
