using .GraphSignal, .GraphPartition, .BasisSpecification, LinearAlgebra

"""
    dmatrix = dvec2dmatrix(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec)

Given a vector of expansion coefficients, convert it to a matrix.

### Input Arguments
* `dvec::Matrix{Float64}`: a vector of expansion coefficients
* `GP::GraphPart`: a GraphPart object
* `BS::BasisSpec`: a BasisSpec object

### Output Argument
* `dmatrix::Array{Float64,3}`: a set of matrices of expansion coefficients
"""
function dvec2dmatrix(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec)

    #
    # 0. Preliminaries
    #
    # extract data
    # MATLAB: [levlist,levlengths] = ExtractData(BS,GP);
    #if isempty(BS.levlengths) && GP != nothing
    #    levlist2levlengths!(GP, BS)
    #end

    levlist = BS.levlist
    #levlengths = BS.levlengths

    # constants
    (N, jmax) = Base.size(GP.rs)
    N = N - 1
    fcols = Base.size(dvec, 2)

    # allocate space
    dmatrix = zeros(N, jmax, fcols)

    #
    # 1. Put the entries in the vector into the correct places in the matrix
    #
    # n = 1
    # for row = 1:length(levlist)
    #     dmatrix[n:(n + levlengths[row] - 1), levlist[row], :] =
    #         dvec[n:(n + levlengths[row] - 1), :]
    #     n += levlengths[row]
    # end
    for i = 1:length(levlist)
        dmatrix[levlist[i][1], levlist[i][2], :] = dvec[i,:]
    end
    #
    # 2. Return the `dmatrix` array
    #
    return dmatrix
end # of function dvec2dmatrix


"""
    (dvec, BS) = dmatrix2dvec(dmatrix::Array{Float64,3}, GP::GraphPart)

Given a matrix of expansion coefficients, convert it to a vector.

### Input Arguments
* `dmatrix::Array{Float64,3}`: matrices of expansion coefficients
* `GP::GraphPart`: an input GraphPart object
* `BS::BasisSpec`: an input BasisSpec object

### Outputs Arguments
* `dvec::Matrix{Float64}`: a vector of expansion coefficients
* `BS::BasisSpec`: an output BasisSpec object
"""
function dmatrix2dvec(dmatrix::Array{Float64,3}, GP::GraphPart)

    # constants
    (N, jmax, fcols) = Base.size(dmatrix)

    # THE BASIS IS NOT SPECIFIED ==> retain nonzero coeffs and specify the basis
    ## 0. Preliminaries
    # allocate/initialize
    dvec = dmatrix[:, jmax, :]
    levlist = jmax * ones(Int, N)

    ## 1. Make a vector of the nonzero basis entries
    for j = (jmax - 1):-1:1
        regioncount = countnz(GP.rs[:, j]) - 1
        for r = 1:regioncount
            indr = GP.rs[r, j]:(GP.rs[r + 1, j] - 1)
            if countnz(dmatrix[indr, j, :]) > 0
                dvec[indr, :] = dmatrix[indr, j, :]
                levlist[GP.rs[r, j]:(GP.rs[r + 1, j] - 1)] .= j
            end
        end
    end
    ## 2. Specify the corresponding basis and return things
    levlist = collect(enumerate(levlist))
    BS = BasisSpec(levlist)
    #levlist2levlengths!(GP, BS)
    return dvec, BS
end # of function dmatrix2dvec (with no BS input)


"""
    dvec = dmatrix2dvec(dmatrix::Array{Float64,3}, GP::GraphPart, BS::BasisSpec)

Given a matrix of expansion coefficients, convert it to a vector.
This function assumes that the input coefficient array `dmatrix` is in the coarse-to-fine format. If `BS.c2f == false`, then this function internally converts `dmatrix` into the fine-to-coarse format. Hence, if one supplies the f2c `dmatrix`, the results become wrong, and the subsequent procedure may result in error.

### Input Arguments
* `dmatrix::Array{Float64,3}`: matrices of expansion coefficients
* `GP::GraphPart`: an input GraphPart object
* `BS::BasisSpec`: an input BasisSpec object

### Outputs Arguments
* `dvec::Matrix{Float64}`: a vector of expansion coefficients
"""
function dmatrix2dvec(dmatrix::Array{Float64,3}, GP::GraphPart, BS::BasisSpec)

    # constants
    (N, jmax, fcols) = Base.size(dmatrix)

    # THE BASIS IS SPECIFIED ==> select the corresponding coefficients
    ## 0. Preliminaries
    # allocate space
    dvec = zeros(N, fcols)
    # if isempty(BS.levlengths)
    #     levlist2levlengths!(GP, BS)
    # else
    #     @warn("We assume that levlengths of the input BS were already computed properly prior to dmatrix2dvec.")
    # end
    levlist = BS.levlist
    #levlengths = BS.levlengths
    BSc2f = BS.c2f
    # put dmatrix in the fine-to-coarse arrangement, if necessary
    if !BSc2f # This assumes that dmatrix has not yet been arranged in f2c
    #if !BSc2f && isempty(GP.rsf2c)
        dmatrix = fine2coarse!(GP, dmatrix = dmatrix, coefp = true)
    end

    ## 1. Make a vector of the matrix entries specified by the BS object
    # n = 1
    # for row = 1:length(levlist)
    #     dvec[n:(n + levlengths[row] - 1), :] =
    #         dmatrix[n:(n + levlengths[row] - 1), levlist[row], :]
    #     n += levlengths[row]
    # end

    for i = 1:length(levlist)
        dvec[i,:] = dmatrix[levlist[i][1], levlist[i][2], :]
    end

    ## 2. Return the dvec
    return dvec
end # of dmatrix2dvec (with BS input)


"""
    levlist2levlengths!(GP::GraphPart, BS::BasisSpec)

Compute the levlengths info for a BasisSpec object

### Input Arguments
* `GP::GraphPart`: a GraphPart object
* `BS::BasisSpec`: BasisSpec object, without `levlengths`; after this function, the `levlengths` field is filled.
"""
function levlist2levlengths!(GP::GraphPart, BS::BasisSpec)

    # 0. Preliminaries
    levlist = BS.levlist
    rs = GP.rs
    # allocate space
    levlengths = zeros(typeof(rs[1]), length(levlist))

    # 1. Determine the length of each basis block specified by levlist
    if BS.c2f                   # coarse-to-fine case
        n = 0
        for row = 1:length(levlist)
            # find the regionstart(s) of the next region
            # MATLAB: IX = find( GP.rs(:,levlist(row))==n+1, 1, 'last' );
            IX = something(findlast(isequal(n+1), rs[:, levlist[row]]), 0)
            #IX = findlast( rs[:, levlist[row]], n + 1 )
            levlengths[row] = rs[IX+1, levlist[row]] - rs[IX, levlist[row]]
            n += levlengths[row]
        end
    else                        # fine-to-coarse case
        n = 0
        if isempty(GP.rsf2c)
            fine2coarse!(GP)
        end
        rsf2c = GP.rsf2c        # This must be here: GP.rsf2c cannot be null!
        for row = 1:length(levlist)
            IX = something(findlast(isequal(n+1), rsf2c[:, levlist[row]]), 0)
            #IX = findlast( rsf2c[:, levlist[row]], n + 1 )
            levlengths[row] = rsf2c[IX + 1, levlist[row]] - rsf2c[IX, levlist[row]]
            n += levlengths[row]
        end
    end

    # 2. Postprocessing
    # get rid of blocks with length 0
    # MATLAB: levlist(levlengths == 0) = [];
    #         levlengths(levlengths == 0) = [];
    #         BS = BasisSpec(levlist,levlengths,c2f,description);
    levlist = levlist[levlengths .!= 0] # . before != is important here!!
    levlengths = levlengths[levlengths .!= 0]
    BS.levlengths = levlengths

end # of function levlist2levlengths!

"""
    (levlistfull, levlengthsfull, transfull) = bsfull(GP::GraphPart, BS::BasisSpec, trans::Vector{Bool})

Given a BasisSpec object, return the full-length, redundant levlist, levlengths, and trans descriptions.

###  Input Arguments
* `GP::GraphPart`: an input GraphPart object
* `BS::BasisSpec`: an input BasisSpec object
* `trans::Vector{Bool}`: a specification of the transforms used for the HGLET-GHWT hybrid transform (default: null)
* `levlengthsp::Bool`: a flag to return levlengthsfull (default: false)
* `transp::Bool`: a flag to return transfull (default: false)

###  Output Arguments
* `levlistfull::Vector{UInt8}`: the full-length, redundant levels list description
* `levlengthsfull::Vector{UInt8}`: the full-length, redundant levels lengths description
* `transfull::Matrix{Bool}`: the full-length, redundant trans description
"""
function bsfull(GP::GraphPart, BS::BasisSpec;
                trans::Vector{Bool} = Vector{Bool}(undef, 0),
                levlengthsp::Bool = false, transp::Bool = false)

    ## 0. Preliminaries

    # extract data
    if isempty(BS.levlengths)
        levlist2levlengths!(GP, BS)
    end
    levlist = BS.levlist
    levlengths = BS.levlengths

    # allocate space
    N = Base.length(GP.ind)
    levlistfull = zeros(Int, N)
    # Assuming that the maximum value of levlist <= jmax \approx log2(N)
    # can be representable by `UInt8`, i.e., N < \approx 5.79E76, which is
    # reasonable. We cannot handle such a large N at this point.

    if levlengthsp
        levlengthsfull = zeros(Int, N)
    end

    if transp && !isempty(trans)
        cols = Base.size(trans, 2)
        # MATLAB: transfull = false(N, cols)
        transfull = falses(N, cols)
    else
        transfull = falses(0)   # transfull is empty.
    end

    ## 1. Fill out the redundant descriptions
    global idx = 0
    for row = 1:length(levlist)
        levlistfull[(idx + 1):(idx + levlengths[row])] .= levlist[row]
        if levlengthsp
            levlengthsfull[(idx + 1):(idx + levlengths[row])] .= levlengths[row]
        end
        if !isempty(transfull)
            transfull[(idx + 1):(idx + levlengths[row]), :] = repmat(trans[row, :], levlengths[row], 1)
        end
        global idx += levlengths[row]
    end

    ## 2. Prepare the returns
    if levlengthsp
        if transp
            return levlistfull, levlengthsfull, transfull
        else
            return levlistfull, levlengthsfull
        end
    else
        if transp
            return levlistfull, transfull
        else
            return levlistfull
        end
    end
end # of function bsfull

"""
    BS = bs_haar(GP::GraphPart)

Specify the Haar basis for a given graph partitioning

### Input Argument
* `GP::GraphPart`: an input GraphPart object

### Output Argument
* `BS::BasisSpec`: a BasisSpec object corresponding to the Haar basis
"""
# function bs_haar(GP::GraphPart)
#
#     # determine jmax
#     jmax = Base.size(GP.rs, 2)
#
#     # allocate space for levlist
#     levlist = zeros(Int, jmax)
#
#     # fill in levlist for the Haar basis
#     levlist[1] = jmax
#     for row = 2:jmax
#         levlist[row] = jmax + 2 - row
#     end
#
#     # make a BasisSpec object
#     BS = BasisSpec(levlist, c2f = false, description = "Haar basis")
#     # fill in the levlengths field of the BasisSpec object
#     levlist2levlengths!(GP, BS)
#     # return it
#     return BS
# end # of function bs_haar

function bs_haar(GP::GraphPart)
    if isempty(GP.tagf2c)
        fine2coarse!(GP)
    end
    # determine jmax
    jmax = Base.size(GP.rs, 2)

    # allocate space for levlist
    levlist = Vector{Tuple{Int,Int}}((findall(x->x==1, GP.tagf2c)))
    pushfirst!(levlist, (1, jmax))

    # make a BasisSpec object
    BS = BasisSpec(levlist, c2f = false, description = "Haar basis")
    # fill in the levlengths field of the BasisSpec object
    #levlist2levlengths!(GP, BS)
    # return it
    return BS
end # of function bs_haar

"""
    BS = bs_level(GP::GraphPart, j::Int, c2f::Bool = true)

Specify the basis corresponding to level j for a given graph partitioning

### Input Arguments
* `GP::GraphPart`: an input GraphPart object
* `j::Int`: the level to which the basis corresponds (`j = 0` is the global level)
* `c2f::Bool`: a flag for c2f or f2c (default: true, i.e., c2f)

### Output Argument
* `BS::BasisSpec`: a BasisSpec object corresponding to the level `j` basis
"""
function bs_level(GP::GraphPart, j::Int; c2f::Bool = true)

    # coarse-to-fine dictionary
    if c2f
        #rspointer = GP.rs[:, j + 1]
        bspec = "coarse-to-fine level $(j)"
    # fine-to-coarse dictionary
    else
        if isempty(GP.rsf2c)
            # if isempty(GP.tag) || isempty(GP.compinfo) # I don't think this part is necessary since GP should be partially filled.
            #    ghwt_core!(GP)
            # end
            fine2coarse!(GP)
        end
        #rspointer = GP.rsf2c[:, j + 1]
        bspec = "fine-to-coarse level $(j)"
    end

    # specify the level j basis
    #Nj = countnz(rspointer) - 1
    if isempty(GP.tag)          # otherwise, it doesn't work for HGLET!
        N = size(GP.rs)[1]-1
    else
        N = size(GP.tag, 1)
    end
    levlist = collect(enumerate((j + 1) * ones(Int,N)))
    BS = BasisSpec(levlist, c2f = c2f, description = bspec)

    # fill in the levlengths field of the BasisSpec object
    #levlist2levlengths!(GP, BS)

    # return it
    return BS
end # of bs_level



"""
    BS = bs_walsh(GP::GraphPart)

Specify the walsh basis corresponding to for a given graph partitioning

### Input Arguments
* `GP::GraphPart`: an input GraphPart object

### Output Argument
* `BS::BasisSpec`: a BasisSpec object corresponding to the walsh basis
"""
function bs_walsh(GP::GraphPart)
    BS = bs_level(GP, 0, c2f = true)
    BS.description = "Walsh basis"
    return BS
end # of bs_level


"""
function dvec_Threshold(dvec::Array{Float64}, SORH::String, keep::Float64,
                GP::GraphPart, BS::BasisSpec)

 Threshold HGLET / GHWT coefficients

### Input Arguments
*   `dvec::Array{Float64,2}`        the vector of expansion coefficients
*   `SORH::String`        use soft ('s') or hard ('h') thresholding
*   `keep::Float64`        a fraction between 0 and 1 which says how many coefficients should be kept
*   `GP::GraphPart`          a GraphPart object, used to identify scaling coefficients
*   `BS::BasisSpec`          a BasisSpec object, used to identify scaling coefficients

### Output Argument
*   `dvec_new::Vector{Float64}`        the thresholded expansion coefficients
*   `kept::Float64`        the number of coefficients kept
"""
function dvec_Threshold(dvec::Array{Float64,2}, SORH::String, keep::Float64,
                GP::GraphPart, BS::BasisSpec)
#
# Apply to 1-D only for now. Need to modify in the future
#
    dvec_new = deepcopy(dvec)
    # Obtain the number of coefficients kept
    if keep > 1 || keep < 0
      error("The input keep should between 0~1")
    else
      kept = Int64(round(keep*size(dvec,1)))
    end

    if SORH == "s" || SORH == "soft"
      ind = sortperm(abs.(dvec[:]), rev = true)
      T = abs(dvec[ind[kept+1]])

      if !isnothing(BS)
        #convert the tag matrix to a vector corresponding to BS
        tag = dmatrix2dvec(Array{Float64,3}(reshape(GP.tag,(size(GP.tag,1),size(GP.tag,2),1))),GP,BS)

        #soft-threshold the non-scaling coefficients
        ind = 1:size(dvec,1)
        ind = ind[(tag .> 0)[:]]
        dvec_new[ind] = sign.(dvec[ind]).*(abs.(dvec[ind]).-T).*(abs.(dvec[ind]).-T .>0)

      # if BS is not given, soft-threshold all coefficients
      else
        dvec_new = sign.(dvec).*(abs.(dvec).-T).*(abs.(dvec).-T.>0)
      end

      kept = count(!iszero, dvec_new)

    elseif SORH == "h" || SORH == "hard"
      ind = sortperm(abs.(dvec[:]), rev = true)
      dvec_new[ind[kept+1:end]] .= 0

      kept = count(!iszero, dvec_new)

    else

      kept = count(!iszero, dvec_new)

    end
    return dvec_new, kept
end
