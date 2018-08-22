using LinearAlgebra
"""
    T = ind_class(N::Int)

    Given `N`, determine what class `ind` (and `rs`) should be

### Input Argument
* `N::Int`: the length of the graph signal

### Output Argument
* `T <: Unsigned`: the unsigned integer class that `ind` should be


Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and revised by Naoki Saito, Feb. 8, 2017
"""
function ind_class(N::Int)
    if N < 2^8-2
        T = UInt8
    elseif N < 2^16-2
        T = UInt16
    elseif N < 2^32-2
        T = UInt32
    elseif N < 2^64-2
        T = UInt64
    else
        T = UInt128             # Hope this does not happen!
    end
    return T
end

"""
    T = tag_class(jmax::Int)

Given `jmax`, determine what class `tag` should be.

### Input Argument
* `jmax::Int`: the number of levels in the recursive partitioning (j = 1:jmax)

### Output Argument
* `T <: Unsigned`: the unsigned integer class that `tag` should be


Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and revised by Naoki Saito, Feb. 8, 2017
"""
function tag_class(jmax::Int)
    if jmax - 1 <= 8
        T = UInt8
    elseif jmax - 1 <= 16
        T = UInt16
    elseif jmax - 1 <= 32
        T = UInt32
    elseif jmax - 1 <= 64
        T = UInt64
    else
        T = UInt128             # Hope this does not happen!
    end
    return T
end

"""
    dmatrix_flatten(dmatrix::Array{Float64,3}, flatten::Any)

Flatten dmatrix using the method specified by the string "flatten"

### Input Arguments
* `dmatrix::Array{Float64,3}`: the matrix of expansion coefficients; after this function is called, it becomes the size of (~, ~, 1).
* `flatten::Any`: the method for flattening dmatrix (see the code in details)
"""
function dmatrix_flatten(dmatrix::Array{Float64,3}, flatten::Any)

    # p-norms or p-quasinorms
    if isa(flatten, Number)
        if flatten == 1
            dmatrix = sum(abs.(dmatrix), dims = 3)
        elseif flatten == 0
            dmatrix = sum(dmatrix != 0, dims = 3)
        else
            dmatrix = (sum(abs.(dmatrix).^flatten, dims = 3)).^(1 / flatten)
        end

        # histogram
        # # # elseif strcmpi(flatten,'histogram') # not yet implemented
        # What I meant here was to implement ASH-like pdf estimation
        # followed by the entropy computation.

        # sum of absolute values = 1-norm
    elseif isa(flatten, Symbol)
        if flatten == :abs
            dmatrix = sum(abs.(dmatrix), 3)
            # standard deviation
        elseif flatten == :std
            # MATLAB: dmatrix = std(dmatrix, 0, 3)
            dmatrix = std(dmatrix, 3)
            # variance
        elseif flatten == :var
            dmatrix = var(dmatrix, 3)
            # inverse variance
        elseif flatten == :invvar
            dmatrix = var(dmatrix, 3)
            dmatrix[ abs.(dmatrix) < eps() ] = eps()
            dmatrix = dmatrix.^-1;
            # max coefficient - min coefficient
        elseif flatten == :minmax || flatten == :maxmin
            # MATLAB: dmatrix = max(dmatrix, [], 3) - min(dmatrix, [], 3)
            dmatrix = maximum(dmatrix, 3) - minimum(dmatrix, 3)
            # sum (no absolute value)
        elseif flatten == :sum
            dmatrix = sum(dmatrix, 3)
            # the sum of the differences of consecutive coeffs
        elseif flatten == :sumdiff
            dmatrix = sum(dmatrix[:, :, 2:end] - dmatrix[:, :, 1:(end - 1)], 3)
            # the sum of the absolute values of the differences of consecutive coeffs
        elseif flatten == :sumabsdiff
            dmatrix = sum(abs.(dmatrix[:, :, 2:end] - dmatrix[:, :, 1:(end - 1)]), 3)
            # Shannon entropy
        elseif flatten === :entropy
            # MATLAB: p = abs(dmatrix)/norm(dmatrix[:], 'fro')
            # p = abs(dmatrix)/norm(dmatrix[:])
            # dmatrix = sum(p. * log2(p), 3)
            # Normally, we flatten to get an energy matrix, and then take entropy
            dmatrix = sum(abs(dmatrix).^2, 3)
            # now, all the column vectors have the same norm, so normalize it.
            dmatrix /= sum(dmatrix[:, 1, 1])
            # now, each column has a unit norm.
            dmatrix = sum(dmatrix .* log2(dmatrix), 3)
            # threshold coefficients below 0.5*norm(dmatrix)/length(dmatrix) and sum
        elseif flatten == :sub
            # MATLAB: t = 0.5 * norm(dmatrix[:], 'fro') / numel(dmatrix)
            t = 0.5 * norm(dmatrix[:]) / length(dmatrix)
            dmatrix[abs.(dmatrix) .< t] = 0
            dmatrix = sum(dmatrix, 3)
            # default (1-norm)
        else
            @warn("the specified flatten symbol $(flatten) is not recognized; hence we assume it as 1-norm.")
            dmatrix = sum(abs.(dmatrix), 3)
        end
    else
        @warn("the specified flatten argument $(flatten) is neither of number nor of symbol type; hence we assume it as 1-norm.")
        dmatrix = sum(abs.(dmatrix), 3)
    end
    return dmatrix
end # of function dmatrix_flatten!

"""
    costfun = cost_functional(cfspec::Any)

Determine the cost functional to be used by the best-basis algorithm.

### Input Argument
* `cfspec::Any`: the specification for the cost functional

### Output Argument
* `costfun::Function`: the cost functional (as a function_handle)
"""
function cost_functional(cfspec::Any)
    if isa(cfspec, Number)
        return function (x) return norm(x, cfspec) end
    elseif isa(cfspec, Function)
        return function (x) return cfspec(x) end
    else return function (x) return norm(x, 1) end # by default 1-norm
        # else return function (x) return norm(x, 0.1) end # by default 0.1-quasinorm
    end
end # of function cost_functional

"""
    (dvec, levlist) = bbchange(dvec::Vector{Float64}, j::Int)

Change to the new best basis

### Input Arguments
* `dvec::Vector{Float64}`: an input coefficient matrix
* `j::Int`: a level index

### Output Arguments
* `dvec::Vector{Float64}`: an output coefficient matrix
* `levlist::Vector{Int}`: an output levlist
"""
function bbchange(dvec::Vector{Float64}, j::Int)
    n = Base.size(dvec, 1) # assume each column corresponding to one signal
    levlist = zeros(Int, n)
    levlist[1] = j
    return dvec, levlist        # might need deepcopy here.
end


"""
    function rs_to_region(rs::Matrix{Any}, tag::Matrix{Any})

From the imformation of rs, tag from GP, compute the tag_r matrix, which have the same
size of dmatrix (expansion coefficient matrix). Each element indicates the place of the
coefficient in the expansion tree.

### Input Arguments
* `rs::Matrix{Any}`: rs from GP, showing information of the partition tree
* `tag::Matrix{Any}`: tag from GP, indicating coefficients tag

### Output Arguments
* `tag_r::Matrix{UInt32}`: showing information of the partition tree, same size as dmatrix
"""
function rs_to_region(rs::Matrix{<:Any}, tag::Matrix{<:Any})
  (m,n) = size(tag)
  tag_r = zeros(m,n)
  for j = 1:(n-1)
    regioncount = count(!iszero, rs[:,j]) - 1
    for r = 1:regioncount
      rs1 = rs[r,j]
      rs3 = rs[r+1,j]
      s = rs3-rs1
      if s==1 #the region and subregion have only one element
        tag_r[rs1,j+1] = 2*tag_r[rs1,j]
      elseif s > 1
        # rs2 marks the start of the second subregion
        rs2 = rs1 + 1
        while rs2 < rs3 && tag[rs2, j+1]!=0
          rs2 = rs2 + 1;
        end

        # the parent region is a copy of the subregion
        if rs2 == rs3
          tag_r[rs1:rs3-1,j+1] = 2*tag_r[rs1:rs3-1,j]

        # the parent region has 2 child regions
        else
        tag_r[rs1:rs2-1,j+1] = 2*tag_r[rs1:rs2-1,j]
        tag_r[rs2:rs3-1,j+1] = 2*tag_r[rs2:rs3-1,j] .+ 1
        end
      end
    end
  end
return Matrix{UInt32}(tag_r)
end


"""
function nonorth2relerror(nonorth::array{Float64,1},B::array{Float64,2})

Given a vector 'nonorth' of non-orthonormal expansion coefficients and
the matrix 'B' such that B*nonorth is the original signal, return a
vector of relative approximation errors when retaining the 1,2,...,N
largest coefficients in magnitude.

### Input argument
* `nonorth`:    a vector of non-orthonormal expansion coefficients
* `B`:          the matrix whose such that B*nonorth is the original signal

### Output argument
* `relerror`:   a vector of relative approximation errors

"""
function nonorth2relerror(nonorth::Array{Float64,1},B::Array{Float64,2})

    # sort the expansion coefficients
    IX = sortperm(nonorth, by = abs, rev = true)

    # generate a matrix where column j contains the j largest coefficients
    matrix0 = full(UpperTriangular(repmat(nonorth[IX],1,length(IX))))
    matrix = deepcopy(matrix0)
    matrix[IX,:] = matrix0

    # the original signal and a matrix of reconstructions
    f = B*nonorth
    recons = B*matrix

    # the relative errors
    relerror = sum((repmat(f,1,length(IX)) - recons).^2,1)'.^(1/2)./norm(f,2)
    return relerror
end


"""
function orth2relerror(orth)

Given a vector 'orth' of orthonormal expansion coefficients, return a
vector of relative approximation errors when retaining the 1,2,...,N
largest coefficients in magnitude.

### Input argument
* `orth`:    a vector of orthonormal expansion coefficients

### Output argument
* `relerror`:   a vector of relative approximation errors

"""
function orth2relerror(orth::Array{Float64,1})
    # sort the coefficients
    orth = sort(orth.^2, rev = true)

    #compute the relative errors
    relerror = ((abs.(sum(orth) .- cumsum(orth))).^(1/2))/sum(orth).^(1/2)

    return relerror
end
