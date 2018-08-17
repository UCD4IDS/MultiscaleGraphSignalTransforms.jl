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
    coeffdict = Array{Dict{Tuple{Int64,Int64},Float64}}(n)

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
    coeffdict_new = Array{Dict{Tuple{Int64,Int64},Float64}}(jmax-1)
    tag_tf = Array{Dict{Tuple{Int64,Int64},Bool}}(jmax-1)

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
    tag_tf_b_new = Array{Dict{Tuple{Int64,Int64},Bool}}(jmax)
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
* `bestbasis`: the values of the coeffcients stored in vector

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function ghwt_tf_bestbasis_new(dmatrix::Matrix{Float64},GP::GraphPart)
    (m,n) = size(dmatrix)

    # Initialization. Store the expanding coeffcients from matrix into a list of dictionary (inbuilt hashmap in Julia)
    # The entry `coeffdict[j][(k,l)]` corresponds to the coefficient of basis-vector on level j with region k and tag l.
    tag = convert(Array{Int64,2},GP.tag)
    tag_r = convert(Array{Int64,2},rs_to_region(GP.rs, GP.tag))
    (m,n) = size(dmatrix)
    coeffdict = Array{Dict{Tuple{Int64,Int64},Float64}}(n)
    for i = 1:n
        coeffdict[i]= Dict{Tuple{Int64,Int64},Float64}((tag_r[j,i],tag[j,i]) => dmatrix[j,i] for j=1:m)
    end

    # TAG_tf stores the time-or-frequency information on every iteration
    # COEFFDICT stores the corresponding cost-functional values.
    COEFFDICT = Vector(n)
    TAG_tf = Vector(n)

    # Initialization of the first iteration
    COEFFDICT[1] = coeffdict
    TAG_tf_init = Array{Dict{Tuple{Int64,Int64},Bool}}(n)
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
    return bestbasis_tag_matrix, bestbasis
end
