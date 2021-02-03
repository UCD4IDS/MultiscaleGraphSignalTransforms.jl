module GHWT

include("utils.jl")

using ..GraphSignal, ..GraphPartition, ..BasisSpecification, LinearAlgebra

include("common.jl")

export ghwt_core!, ghwt_analysis!, fine2coarse!, ghwt_synthesis, ghwt_c2f_bestbasis, ghwt_f2c_bestbasis, ghwt_bestbasis, GHWT_jkl

"""
    ghwt_core!(GP::GraphPart)

performs the forward GHWT transform, which computes tag and compinfo of GP, but not the expansion coefficients.

### Input Arguments
* `GP::GraphPart`: a GraphPart object (with or without tag and compinfo data); note that tag and compinfo will be modified after this function is executed
"""
function ghwt_core!(GP::GraphPart)

    #
    # 0. Preliminaries
    #
    rs = GP.rs                  # for simple notation!
    (N, jmax) = Base.size(rs)
    N = N - 1

    # tag -- when expressed in binary, tag indicates the sequence of
    # average (0's) and difference (1's) operations used to obtain the given
    # basis vector
    if isempty(GP.tag)
        GP.tag = zeros(Int, N, jmax)
    end
    tag = GP.tag                # for simple notation!

    # allocate space for compinfo (used for synthesis from f2c dictionary)
    # if tag == 0 && n == 1:                        compinfo = 0
    # if tag == 0 && n >= 2:                        compinfo = n1
    # if tag == 1:                                  compinfo = n2
    # if tag >= 2 && coeff is formed from 2 coeffs: compinfo = 1
    # if tag >= 2 && coeff is formed from 1 coeff:  compinfo = 0
    if isempty(GP.compinfo)
        GP.compinfo = zeros(Int, N, jmax)
    end
    compinfo = GP.compinfo      # for simple notation!

    #
    # 1. Perform the transform
    #
    for j = (jmax - 1):-1:1 # # Bottom-1 (jmax-1) to up to the coarsest (j=1)
        regioncount = count(!iszero,rs[:, j]) - 1
        for r = 1:regioncount
            rs1 = rs[r, j]      # the start of the 1st subregion
            rs3 = rs[r + 1, j]  # 1 + the end of the 2nd subregion
            n = rs3 - rs1       # # of points in the current region
            if n > 1            # single node region passes this part
                rs2 = rs1 + 1   # the start of the 2nd subregion
                while rs2 < rs3 && tag[rs2, j + 1] != 0 ### && rs2 < N+1
                    rs2 += 1
                end
                if rs2 == rs3 # the parent region is a copy of the subregion
                    tag[rs1:(rs3 - 1), j] = tag[rs1:(rs3 - 1), j + 1]
                else           # the parent region has 2 child regions
                    n1 = rs2 - rs1 # # of pts in the 1st subregion
                    n2 = rs3 - rs2 # # of pts in the 2nd subregion
                    # SCALING COEFFICIENT (n > 1)
                    compinfo[rs1, j] = n1;
                    # HAAR COEFFICIENT
                    compinfo[rs1 + 1, j] = n2; tag[rs1 + 1, j] = 1
                    # WALSH COEFFICIENTS
                    # sweep through the coefficients in subregions 1 & 2
                    parent = rs1 + 2 # ind of the new coeffs to be created on level j
                    child1 = rs1 + 1 # ind of the current coeffs in subregion 1
                    child2 = rs2 + 1 # ind of the current coeffs in subregion 2
                    while parent < rs3
                        # no matching coefficient (use subregion 1)
                        if child1 < rs2 && (child2 == rs3 || tag[child1, j + 1] < tag[child2, j + 1])
                            tag[parent, j] = 2 * tag[child1, j + 1]
                            child1 += 1; parent += 1
                        # no matching coefficient (use subregion 2)
                        elseif child2 < rs3 && (child1 == rs2 || tag[child2, j + 1] < tag[child1, j + 1])
                            tag[parent, j] = 2 * tag[child2, j + 1]
                            child2 += 1; parent += 1
                        # matching coefficients
                        else
                            tag[parent, j] = 2 * tag[child1, j + 1]
                            #tag[parent + 1, j] = 2 * tag[child1, j + 1] + 1
                            tag[parent + 1, j] = tag[parent, j] + 1
                            compinfo[parent, j] = 1; compinfo[parent + 1, j]  = 1
                            child1 += 1; child2 += 1; parent += 2
                        end # of if child1 < rs2 ...
                    end # of while
                end # of if rs2 == rs3 ...
            elseif n <= 0 # n must be positive
                error("n must be positive: n = ", n)
            end # of if n > 1 ... elseif n <= 0 ...; note: do nothing for n==1
        end # of for r = 1:regioncount
    end # of for j = (jmax - 1):-1:1
end # of function ghwt_core!


"""
    ghwt_core!(GP::GraphPart, dmatrix::Array{Float64,3})

performs the forward GHWT transform, which generates expansion coefficients, tag, and compinfo.

### Input Arguments
* `GP::GraphPart`: a GraphPart object (with or without tag and compinfo data); note that tag and compinfo will be modified after this function is executed
* `dmatrix::Array{Float64,3}`: a matrix of expansion coefficients with only the last column filled in as the original input signal(s) `f`; after this function is executed, this matrix contains whole expansion coefficients of the input signal(s) relative to the GHWT dictionary
"""
function ghwt_core!(GP::GraphPart, dmatrix::Array{Float64,3})

    #
    # 0. Preliminaries
    #
    rs = GP.rs                  # for simple notation!
    (N, jmax) = Base.size(rs)
    N = N - 1

    # tag -- when expressed in binary, tag indicates the sequence of
    # average (0's) and difference (1's) operations used to obtain the given
    # basis vector
    if isempty(GP.tag)
        GP.tag = zeros(Int, N, jmax)
    end
    tag = GP.tag                # for simple notation!

    # allocate space for compinfo (used for synthesis from fine-to-coarse dictionary)
    # if tag == 0 && n == 1:                        compinfo = 0
    # if tag == 0 && n >= 2:                        compinfo = n1
    # if tag == 1:                                  compinfo = n2
    # if tag >= 2 && coeff is formed from 2 coeffs: compinfo = 1
    # if tag >= 2 && coeff is formed from 1 coeff:  compinfo = 0
    if isempty(GP.compinfo)
        GP.compinfo = zeros(Int, N, jmax)
    end
    compinfo = GP.compinfo      # for simple notation!

    #
    # 1. Perform the transform
    #
    for j = (jmax - 1):-1:1 # Bottom-1 (jmax-1) to up to the coarsest (j=1)
        regioncount = count(!iszero,rs[:, j]) - 1
        for r = 1:regioncount
            rs1 = rs[r, j]      # the start of the 1st subregion
            rs3 = rs[r + 1, j]  # 1 + the end of the 2nd subregion
            n = rs3 - rs1       # # of points in the current region
            if n == 1           # single node region
                # SCALING COEFFICIENT (n == 1)
                dmatrix[rs1, j, :] = dmatrix[rs1, j + 1, :]
            elseif n > 1
                rs2 = rs1 + 1 # the start of the 2nd subregion
                while rs2 < rs3 && tag[rs2, j + 1] != 0 ### && rs2 < N+1
                    rs2 += 1
                end
                if rs2 == rs3 # the parent region is a copy of the subregion
                    dmatrix[rs1:(rs3 - 1), j, :] = dmatrix[rs1:(rs3 - 1), j + 1, :]
                    tag[rs1:(rs3 - 1), j] = tag[rs1:(rs3 - 1), j + 1]
                else       # the parent region has 2 child regions
                    n1 = rs2 - rs1 # # of pts in the 1st subregion
                    n2 = rs3 - rs2 # # of pts in the 2nd subregion

                    # SCALING COEFFICIENT (n > 1)
                    dmatrix[rs1, j, :] =
                        ( sqrt(n1) * dmatrix[rs1, j + 1, :] +
                          sqrt(n2) * dmatrix[rs2, j + 1, :] ) / sqrt(n)
                    compinfo[rs1, j] = n1

                    # HAAR COEFFICIENT
                    dmatrix[rs1 + 1, j, :] =
                        ( sqrt(n2) * dmatrix[rs1, j + 1, :] -
                          sqrt(n1) * dmatrix[rs2, j + 1, :] ) / sqrt(n)
                    # Note these sqrt's arguments are not typos!!
                    compinfo[rs1 + 1, j] = n2; tag[rs1 + 1, j] = 1

                    # WALSH COEFFICIENTS
                    # sweep through the coefficients in subregions 1 & 2
                    parent = rs1 + 2 # ind of the new coeffs to be created on level j
                    child1 = rs1 + 1 # ind of the current coeffs in subregion 1
                    child2 = rs2 + 1 # ind of the current coeffs in subregion 2
                    while parent < rs3
                        # no matching coefficient (use subregion 1)
                        if child1 < rs2 && (child2 == rs3 || tag[child1, j + 1] < tag[child2, j + 1])
                            dmatrix[parent, j, :] = dmatrix[child1, j + 1, :]
                            tag[parent, j] = 2 * tag[child1, j + 1]
                            child1 += 1; parent += 1
                        # no matching coefficient (use subregion 2)
                        elseif child2 < rs3 && (child1 == rs2 || tag[child2, j + 1] < tag[child1, j + 1])
                            dmatrix[parent, j, :] = dmatrix[child2, j + 1, :]
                            tag[parent, j] = 2 * tag[child2, j + 1]
                            child2 += 1; parent += 1
                        # matching coefficients
                        else
                            dmatrix[parent, j, :] =
                                ( dmatrix[child1, j + 1, :] +
                                  dmatrix[child2, j + 1, :] ) / sqrt(2)
                            tag[parent, j] = 2 * tag[child1, j + 1]
                            compinfo[parent, j] = 1

                            dmatrix[parent + 1, j, :] =
                                ( dmatrix[child1, j + 1, :] -
                                  dmatrix[child2, j + 1, :] ) / sqrt(2)
                            #tag[parent + 1, j] = 2 * tag[child1, j + 1] + 1
                            tag[parent + 1, j] = tag[parent, j] + 1
                            compinfo[parent + 1, j]  = 1
                            child1 += 1; child2 += 1; parent += 2
                        end # of if child1 < rs2 ...
                    end # of while
                end # of if rs2 == rs3 ...
            else # n must be positive
                error("n must be positive: n = ", n)
            end # of if n == 1 ... elseif n > 1 ...
        end # of for r = 1:regioncount
    end # of for j = (jmax - 1):-1:1
end # of function ghwt_core!


"""
    dmatrix = ghwt_analysis!(G::GraphSig, GP::GraphPart = nothing, c2f::Bool = true)

For a GraphSig object `G`, generate the matrix of GHWT expansion coefficients

### Input Arguments
* `G::GraphSig`: an input GraphSig object
* `GP::GraphPart`: an input GraphPart object (optional); after this function is run, `GP`'s `compinfo`, `tag`, etc. are filled

### Output Argument
* `dmatrix::Array{Float64,3}`: the 3D array of expansion coefficients (i.e., for each input signal vector, the matrix of coefficients; hence, for multiple input signals, the coefficients are organized as a 3D array)
"""
function ghwt_analysis!(G::GraphSig; GP::GraphPart = nothing, c2f::Bool = true)

    #
    # 0. Preliminaries
    #
    if GP == nothing
        GP = partition_tree_fiedler(G)
    end

    # Simplify the subsequent notation without overhead. We can do this
    # since `ind`, `rs`, `f` share the same memory as `GP.ind`, `GP.rs`, `G.f`
    ind = GP.ind
    rs = GP.rs
    f = G.f
    (N, jmax) = Base.size(rs)
    N = N - 1

    # allocate space for the expansion coefficients
    fcols = size(f, 2)
    dmatrix = zeros(N, jmax, fcols)
    dmatrix[:, jmax, :] = f[ind, :]

    ## 1. Perform the transform

    # generate expansion coefficients, tag, and compinfo
    ghwt_core!(GP, dmatrix)

    # fill in the remaining (i.e. fine-to-coarse) GraphPart fields
    if c2f
        return dmatrix
    else
        return fine2coarse!(GP, dmatrix=dmatrix, coefp = true)
    end

end # of ghwt_analysis!


"""
    (dmatrixf2c, IX) = fine2coarse!(GP::GraphPart;
    dmatrix::Array{Float64,3} = zeros(0, 0, 0),
    coefp::Bool = false, indp::Bool = false)

Fill in the fine-to-coarse info (rs2f2c, tagf2c, and compinfof2c) in a
GraphPart object.  Also, rearrange a matrix of expansion coefficients.

### Input Arguments
* `GP::GraphPart`: an input GraphPart object without fine-to-coarse info (rsf2c, tagf2c, compinfof2c); after this function, rsf2c, tagf2c, compinfof2c are filled. Note that rs, tag, compinfo are intact.
* `dmatrix::Array{Float64,3}`: a matrix of expansion coefficients in coarse-to-fine arrangement (default: null matrix)
* `coefp::Bool`: a flag to return the rearranged f2c coefficients (default: false)
* `indp::Bool`: a flag to return the reordering index (default: false)

### Output Arguments
* `dmatrixf2c::Array{Float64,3}`: a matrix of expansion coefficients in fine-to-coarse arrangement (if requested)
* `IX::Vector{Unsigned}`: the reordering index for all levels
"""
function fine2coarse!(GP::GraphPart;
                      dmatrix::Array{Float64,3} = zeros(0, 0, 0),
                      coefp::Bool = false, indp::Bool = false)
    #
    # 0. Preliminaries
    #

    # get constants
    (N, jmax) = Base.size(GP.rs)
    N = N - 1

    # allocate spaces
    if isempty(GP.rsf2c)
        GP.rsf2c = zeros(Int, N + 1, jmax)
    end
    GP.rsf2c[1, :] .= 1
    GP.rsf2c[2, 1] = N + 1
    if isempty(GP.tagf2c)
        GP.tagf2c = zeros(Int, N, jmax)
    end
    if isempty(GP.compinfof2c)
        GP.compinfof2c = zeros(Int, N, jmax)
    end
    IX = zeros(Int, N, jmax)
    if coefp
        dmatrixf2c = zeros(Base.size(dmatrix))
    end

    # make sure the coarse-to-fine tag field is filled in
    if isempty(GP.tag) || isempty(GP.compinfo)
        ghwt_core!(GP)
    end

    # for simple notation!
    ind = GP.ind; rs = GP. rs
    tag = GP.tag; compinfo = GP.compinfo
    rsf2c = GP.rsf2c; tagf2c = GP.tagf2c; compinfof2c = GP.compinfof2c

    #
    # 1. Generate the fine-to-coarse dictionary
    #

    # put the coefficients into fine-to-coarse order
    for j = 1:jmax
        # put the basis into tag order
        # MATLAB: [GP.tagf2c(:,j),IX(:,j)] = sort(GP.tag(:,jmax+1-j));
        IX[:, j] = sortperm(tag[:, jmax + 1 - j])
        tagf2c[:, j] = tag[IX[:, j], jmax + 1 - j]
        compinfof2c[:, j] = compinfo[IX[:, j], jmax + 1 - j]

        if coefp
            dmatrixf2c[:, j, :] = dmatrix[IX[:, j], jmax + 1 - j, :]
        end

        # fill in the fine-to-coarse regionstarts
        if j > 1
            rsf2crow = 2
            for row = 2:N
                if tagf2c[row, j] > tagf2c[row - 1, j]
                    rsf2c[rsf2crow, j] = row
                    rsf2crow += 1
                end
            end
            rsf2c[rsf2crow, j] = N + 1
        end
    end # of for j = 1:jmax

    #
    # 2. Return the appropriate variables if requested
    #
    if coefp && indp
        return dmatrixf2c, IX
    elseif coefp && !indp
        return dmatrixf2c
    elseif !coefp && indp
        return IX
    end
end # of function fine2coarse!


"""
    f = ghwt_synthesis(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec)

Given a vector of GHWT expansion coefficients and info about the graph
partitioning and the choice of basis, reconstruct the signal

### Input Arguments
* `dvec::Matrix{Float64}`: the expansion coefficients corresponding to the chosen basis
* `GP::GraphPart`: an input GraphPart object
* `BS::BasisSpec`: an input  BasisSpec object

### Output Arguments
* `f::Matrix{Float64}`: the reconstructed signal(s)
"""
function ghwt_synthesis(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec)

    #
    # 0. Preliminaries
    #
    # if necessary, fill in the tag info in GP
    if isempty(GP.tag)
        ghwt_core!(GP)
    end

    # figure out which dictionary is used: coarse-to-fine or fine-to-coarse
    if !BS.c2f && isempty(GP.tagf2c)
        fine2coarse!(GP)
    end

    # constants
    N = Base.length(GP.ind)
    jmax = Base.size(GP.rs, 2)
    fcols = Base.size(dvec, 2)

    # fill in the appropriate entries of dmatrix
    dmatrix = dvec2dmatrix(dvec, GP, BS)

    #
    # 1. Synthesis from the basis dictionary
    #
    if BS.c2f                   # for coarse-to-fine case
        # for simpler notation!
        dvec_loc = zeros(Int, size(GP.tag))
        for i = 1:length(BS.levlist)
            dvec_loc[BS.levlist[i][1], BS.levlist[i][2]] = 1
        end

        rs = GP.rs
        tag = GP.tag
        for j = 1:(jmax - 1)    # from top to bottom-1
            regioncount = count(!iszero, rs[:, j]) - 1
            for r = 1:regioncount
                rs1 = rs[r, j]      # the start of the 1st subregion
                rs3 = rs[r + 1, j]  # 1 + the end of the 2nd subregion
                n = rs3 - rs1       # # of points in the current region
                # only proceed forward if coefficients do not exist
                #if count(!iszero, dmatrix[rs1:(rs3 - 1), j + 1, :]) == 0 &&
                #    count(!iszero, dmatrix[rs1:(rs3 - 1), j, :]) > 0
                if count(!iszero, dvec_loc[rs1:(rs3-1),j]) > 0
                    if n == 1   # single node region
                        # SCALING COEFFICIENT (n == 1)
                        if dvec_loc[rs1,j] == 1
                            dmatrix[rs1, j + 1, :] = dmatrix[rs1, j, :]
                            dvec_loc[rs1,j+1] = 1
                        end
                    elseif n > 1
                        rs2 = rs1 + 1 # the start of the 2nd subregion
                        while rs2 < rs3 && tag[rs2, j + 1] != 0 ### && rs2 < N+1
                            rs2 += 1
                        end
                        if rs2 == rs3 # the parent is a copy of the subregion
                            if dvec_loc[rs1:rs3-1,j] == 1
                                dmatrix[rs1:(rs3 - 1), j + 1, :] = dmatrix[rs1:(rs3 - 1), j, :]
                                dvec_loc[rs1:rs3-1,j+1] = 1
                            end
                        else   # the parent region has 2 child regions
                            n1 = rs2 - rs1 # # of pts in the 1st subregion
                            n2 = rs3 - rs2 # # of pts in the 2nd subregion

                            # SCALING COEFFICIENTS (n > 1)
                            if dvec_loc[rs1,j] == 1 && dvec_loc[rs1+1,j] == 1
                                dmatrix[rs1, j + 1, :] =
                                    ( sqrt(n1) * dmatrix[rs1, j, :] +
                                      sqrt(n2) * dmatrix[rs1 + 1, j, :] ) / sqrt(n)

                            # HAAR COEFFICIENTS
                                dmatrix[rs2, j + 1, :] =
                                    ( sqrt(n2) * dmatrix[rs1, j, :] -
                                      sqrt(n1) * dmatrix[rs1 + 1, j, :] ) / sqrt(n)

                                dvec_loc[rs1,j+1] = 1
                                dvec_loc[rs2,j+1] = 1
                            end

                            # WALSH COEFFICIENTS
                            # search through the remaining coefs in each subregion
                            parent = rs1 + 2
                            child1 = rs1 + 1
                            child2 = rs2 + 1
                            while child1 < rs2 || child2 < rs3
                                # subregion 1 has the smaller tag
                                if child2 == rs3 ||
                                    (tag[child1, j + 1] < tag[child2, j + 1] &&
                                     child1 < rs2)

                                    if dvec_loc[parent,j]==1
                                        dmatrix[child1, j + 1, :] =
                                            dmatrix[parent, j, :]
                                        dvec_loc[child1, j+1] =1
                                    end
                                    child1 += 1; parent += 1

                                # subregion 2 has the smaller tag
                                elseif child1 == rs2 ||
                                    (tag[child2, j + 1] < tag[child1, j + 1] &&
                                     child2 < rs3)

                                    if dvec_loc[parent, j] == 1
                                        dmatrix[child2, j + 1, :] =
                                            dmatrix[parent, j, :]
                                        dvec_loc[child2, j+1] = 1
                                    end
                                    child2 += 1; parent += 1

                                # both subregions have the same tag
                                else
                                    if dvec_loc[parent,j] == 1 && dvec_loc[parent+1, j] == 1
                                        dmatrix[child1, j + 1, :] =
                                            ( dmatrix[parent, j, :] +
                                              dmatrix[parent + 1, j, :] ) / sqrt(2)
                                        dmatrix[child2, j + 1, :] =
                                            ( dmatrix[parent, j, :] -
                                              dmatrix[parent + 1, j, :] ) / sqrt(2)
                                        dvec_loc[child1,j+1]=1
                                        dvec_loc[child2,j+1]=1
                                    end
                                    child1 += 1; child2 += 1; parent += 2
                                end # of if child2 == r3 ...
                            end # of while child1 < rs2 ...
                        end # of if rs2 == rs3 ... else
                    end # of if n == 1 ... elseif n > 1 ...
                end # of if countnz(...)
            end # of for r = 1:regioncount
        end # of for j = 1:(jmax - 1)
        # MATLAB: if fcols == 1
        #            f = dmatrix[:, end]
        #         else
        #            f = reshape(dmatrix[:, end, :], N, fcols)
        #         end
        # GS.f = dmatrix[:, end, :] # Do this later for clarity
    else                        # for fine-to-coarse case
        # for simpler notation!
        rsf2c = GP.rsf2c
        tagf2c = GP.tagf2c
        compinfof2c = GP.compinfof2c
        for j = (jmax - 1):-1:1 # from bottom-1 to top
            regioncount = count(!iszero, rsf2c[:, j]) - 1
            for r = 1:regioncount
                rs1 = rsf2c[r, j] # the start of the 1st subregin
                rs3 = rsf2c[r + 1, j] # 1 + the end of the 2nd subregion
                # only proceed forward if coefficients do not exist
                if count(!iszero, dmatrix[rs1:(rs3 - 1), j, :]) == 0 &&
                    count(!iszero, dmatrix[rs1:(rs3 - 1), j + 1, :]) > 0
                    # one subregion ==> copy the coefficients
                    if tagf2c[rs1, j + 1] == tagf2c[rs3 - 1, j + 1]
                        dmatrix[rs1:(rs3 - 1), j, :] =
                            dmatrix[rs1:(rs3 - 1), j + 1, :]
                    # two subregions ==> compute the coefficients
                    else
                        rs2 = rs1 + 1 # the start of the 2nd subregion
                        while rs2 < rs3 &&
                            tagf2c[rs1, j + 1] == tagf2c[rs2, j + 1] ### && rs2 < N+1
                            rs2 += 1
                        end
                        if rs2 == rs3 # the parent is a copy of the subregion
                            # no need to do the following
                            # dmatrix[rs1:(rs3 - 1), (j + 1) , :] =
                            #  dmatrix[rs1:(rs3 - 1) , j, :]
                            nothing
                        else   # the parent region has 2 child regions
                            parent = rs1 # the current coef in the parent
                            child1 = rs1 # the current coefs in the children
                            child2 = rs2

                            # SCALING & HAAR COEFFICIENTS
                            if tagf2c[rs1, j] == 0
                                while parent < rs3
                                    # 1 coef formed from 1 coef
                                    if compinfof2c[child1, j + 1] == 0
                                        dmatrix[parent, j, :] =
                                            dmatrix[child1, j + 1, :]
                                        parent += 1; child1 += 1
                                    # 2 coefs formed from 2 coefs
                                    else
                                        n1 = compinfof2c[child1, j + 1]
                                        n2 = compinfof2c[child2, j + 1]
                                        n = n1 + n2 # # of points in the parent
                                        # SCALING COEFFICIENT
                                        dmatrix[parent, j, :] =
                                            ( sqrt(n1) * dmatrix[child1, j + 1, :] +
                                              sqrt(n2) * dmatrix[child2, j + 1, :] ) / sqrt(n)
                                        # HAAR COEFFICIENT
                                        dmatrix[parent + 1, j, :] =
                                            ( sqrt(n2) * dmatrix[child1, j + 1, :] -
                                              sqrt(n1) * dmatrix[child2, j + 1, :] ) / sqrt(n)
                                        parent += 2; child1 += 1; child2 += 1
                                    end
                                end # of while
                            # WALSH COEFFICIENTS
                            else
                                while parent < rs3
                                    # 1 coef formed from 1 coef
                                    if compinfof2c[child1, j + 1] == 0
                                        dmatrix[parent, j, :] =
                                            dmatrix[child1, j + 1, :]
                                        parent += 1; child1 += 1

                                    # 2 coefs formed from 2 coefs
                                    else
                                        dmatrix[parent, j, :] =
                                            ( dmatrix[child1, j + 1, :] +
                                              dmatrix[child2, j + 1, :] ) / sqrt(2)
                                        dmatrix[parent + 1, j, :] =
                                            ( dmatrix[child1, j + 1, :] -
                                              dmatrix[child2, j + 1, :] ) / sqrt(2)
                                        parent += 2; child1 += 1; child2 += 1
                                    end # of if compinfof2c[...]
                                end # of while parent < rs3
                            end # of if tagf2c[rs1, j] == 0 ... else ...
                        end # of if rs2 == rs3 ... else ...
                    end # of if tagf2c[rs1, j + 1] == tagf2c[rs3 - 1, j + 1] ...
                end # of if countnz(dmatrix[rs1:(rs3 - 1), j, :]) == 0 ...
            end # for r = 1:regioncount
        end # for j=(jmax - 1):-1:1
        # GS.f = dmatrix[:, 1, :]
    end # of if BS.c2f ... else ...

    #
    # 2. Final prepartion for returning variables
    #

    ftemp = dmatrix[:, 1, :] # for fine-to-coarse and make f available
    if BS.c2f
        ftemp = dmatrix[:, end, :] # for coarse-to-fine
    end
    # put the reconstructed values in the correct order and return it.
    f = zeros(size(ftemp))
    f[GP.ind, :] = ftemp

    return f
end # of function ghwt_synthesis (without G input)


"""
    (f, GS) = ghwt_synthesis(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec, G::GraphSig)

Given a vector of GHWT expansion coefficients and info about the graph
partitioning and the choice of basis, reconstruct the signal

### Input Arguments
* `dvec::Matrix{Float64}`: the expansion coefficients corresponding to the chosen basis
* `GP::GraphPart`: an input GraphPart object
* `BS::BasisSpec`: an input  BasisSpec object
* `G::GraphSig`: an input  GraphSig object

### Output Arguments
* `f::Matrix{Float64}`: the reconstructed signal(s)
* `GS::GraphSig`: the reconstructed GraphSig object
"""
function ghwt_synthesis(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec, G::GraphSig)

    # call the normal ghwt_synthesis
    f = ghwt_synthesis(dvec, GP, BS)

    # create a GraphSig object with the reconstructed data and return it with f.
    GS = deepcopy(G)        # generate a new copy of G as GS
    replace_data!(GS, f)    # then replace the graph signal field
    return f, GS
end # of ghwt_synthesis (with G input)

"""
    (dvecc2f, BSc2f) = ghwt_c2f_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart; cfspec::Any = 1.0, flatten::Any = 1.0)

Select the coarse-to-fine best basis from the matrix of GHWT expansion coefficients

### Input Arguments
* `dmatrix::Array{Float64,3}`: the matrix of expansion coefficients
* `GP::GraphPart`: an input GraphPart object
* `cfspec::Any`: the specification of cost functional to be used (default = 1.0, i.e., 1-norm)
* `flatten::Any`: the method for flattening vector-valued data to scalar-valued data (default = 1.0, i.e, 1-norm)

###  Output Arguments
* `dvecc2f::Matrix{Float64}`: the vector of expansion coefficients corresponding to the coarse-to-fine best basis
* `BSc2f::BasisSpec`: a BasisSpec object which specifies the coarse-to-fine best basis
"""
function ghwt_c2f_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart;
                            cfspec::Any = 1.0, flatten::Any = 1.0, j_start::Int = 1, j_end::Int = size(dmatrix,2))

    # determine the cost functional to be used
    costfun = cost_functional(cfspec)

    # constants and dmatrix cleanup
    (N, jmax, fcols) = Base.size(dmatrix)
    if ~(1 <= j_start <= j_end <= jmax)
        error("j_start and j_end should satisfy the inequality 1 <= j_start <= j_end <= "*string(jmax))
    end
    dmatrix[ abs.(dmatrix) .< 10^2 * eps() ] .= 0

    # "flatten" dmatrix
    if fcols > 1
        dmatrix0 = deepcopy(dmatrix)      # keep the original dmatrix as dmatrix0
        if flatten != nothing
            dmatrix = dmatrix_flatten(dmatrix, flatten)
        end
    end

    ## Find the best-basis from the GHWT coarse-to-fine dictionary
    # allocate/initialize
    dvecc2f = dmatrix[:, j_end, 1]
    levlistc2f = j_end * ones(Int, N)

    # set the tolerance
    tol = 10^4 * eps()

    # perform the basis search
    for j = (j_end - 1):-1:j_start
        regioncount = count(!iszero, GP.rs[:, j]) - 1
        for r = 1:regioncount
            indr = GP.rs[r, j]:(GP.rs[r + 1, j]-1)
            ##### compute the cost of the current best basis
            costBB = costfun(dvecc2f[indr])
            costNEW = costfun(dmatrix[indr, j, 1])
            if costBB >= costNEW - tol
                #(dvecc2f[indr], levlistc2f[indr]) = bbchange(dmatrix[indr, j, 1], j)
                dvecc2f[indr] = dmatrix[indr, j, 1]
                levlistc2f[indr] .= j
            end
        end
    end
    #levlistc2f = Array{Int}(levlistc2f[ levlistc2f .!= 0 ])
    levlistc2f = collect(enumerate(levlistc2f))
    BSc2f = BasisSpec(levlistc2f, c2f = true, description = "GHWT c2f Best Basis")
    #levlist2levlengths!(GP, BSc2f)

    # if we flattened dmatrix, then "unflatten" the expansion coefficients
    if fcols > 1
        dvecc2f = dmatrix2dvec(dmatrix0, GP, BSc2f)
    else
        dvecc2f = reshape(dvecc2f,(length(dvecc2f),1))
    end

    # Return things
    return dvecc2f, BSc2f
end # of function ghwt_c2f_bestbasis

"""
    (dvecf2c, BSf2c) = ghwt_f2c_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart; cfspec::Any = 1.0, flatten::Any = 1.0)

Select the fine-to-coarse best basis from the matrix of GHWT expansion coefficients

### Input Arguments
* `dmatrix::Array{Float64,3}`: the matrix of expansion coefficients
* `GP::GraphPart`: an input GraphPart object
* `cfspec::Any`: the specification of cost functional to be used (default: 1.0, i.e., 1-norm)
* `flatten::Any`: the method for flattening vector-valued data to scalar-valued data (default: 1.0, i.e., 1-norm)

###  Output Arguments
* `dvecf2c::Matrix{Float64}`: the vector of expansion coefficients corresponding to the fine-to-coarse best basis
*  `BSf2c::BasisSpec`: a BasisSpec object which specifies the fine-to-coarse best basis
"""
function ghwt_f2c_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart;
                            cfspec::Any = 1.0, flatten::Any = 1.0, j_start::Int = 1, j_end::Int = size(dmatrix,2))

    # determine the cost functional to be used
    costfun = cost_functional(cfspec)

    # constants and dmatrix cleanup
    (N, jmax, fcols) = Base.size(dmatrix)
    if ~(1 <= j_start <= j_end <= jmax)
        error("j_start and j_end should satisfy the inequality 1 <= j_start <= j_end <= "*string(jmax))
    end

    dmatrix[ abs.(dmatrix) .< 10^2 * eps() ] .= 0

    # "flatten" dmatrix
    if fcols > 1
        dmatrix0 = deepcopy(dmatrix)      # keep the original dmatrix as dmatrix0
        if flatten != nothing
            dmatrix = dmatrix_flatten(dmatrix, flatten)
        end
    end

    ## Find the best-basis from the GHWT fine-to-coarse dictionary

    # generate the fine-to-coarse GraphPart fields and coefficient matrix
    dmatrixf2c = fine2coarse!(GP, dmatrix = dmatrix, coefp = true)

    # allocate/initialize
    dvecf2c = dmatrixf2c[:, j_end, 1]
    levlistf2c = j_end * ones(Int, N)

    # set the tolerance
    tol = 10^4 * eps()

    # perform the basis search
    for j = (j_end - 1):-1:j_start
        regioncount = count(!iszero, GP.rsf2c[:, j]) - 1
        for r = 1:regioncount
            indr = GP.rsf2c[r, j]:(GP.rsf2c[r + 1, j] - 1)
            ##### compute the cost of the current best basis
            costBB = costfun(dvecf2c[indr])
            costNEW = costfun(dmatrixf2c[indr, j, 1])
            if costBB >= costNEW - tol
                #(dvecf2c[indr], levlistf2c[indr]) = bbchange(dmatrixf2c[indr, j, 1], j)
                dvecf2c[indr] = dmatrixf2c[indr, j, 1]
                levlistf2c[indr] .= j
            end
        end
    end

    #levlistf2c = Array{Int}(levlistf2c[ levlistf2c .!= 0 ])
    levlistf2c = collect(enumerate(levlistf2c))
    BSf2c = BasisSpec(levlistf2c, c2f = false, description = "GHWT f2c Best Basis")
    #levlist2levlengths!(GP, BSf2c)
    # costf2c = costfun(dvecf2c)

    # if we flattened dmatrix, then "unflatten" the expansion coefficients
    if fcols > 1
        dvecf2c = dmatrix2dvec(dmatrix0, GP, BSf2c)
    else
        dvecf2c = reshape(dvecf2c,(length(dvecf2c),1))
    end

    # Return things
    return dvecf2c, BSf2c
end # of function ghwt_f2c_bestbasis

"""
    (dvec, BS) = ghwt_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart; cfspec::Any, flatten::Any = 1.0)

Select the overall best basis among the c2f and f2c best bases

### Input Arguments
* `dmatrix::Array{Float64,3}`: the matrix of expansion coefficients
* `GP::GraphPart`: an input GraphPart object
* `cfspec::Any`: the specification of cost functional to be used (default: 1.0, i.e., 1-norm)
* `flatten::Any`: the method for flattening vector-valued data to scalar-valued data (default: 1.0, i.e., 1-norm)

###  Output Arguments
* `dvec::Matrix{Float64}`: the vector of expansion coefficients corresponding to the best basis
* `BS::BasisSpec`: a BasisSpec object which specifies the best basis
"""
function ghwt_bestbasis(dmatrix::Array{Float64,3}, GP::GraphPart;
                        cfspec::Any = 1.0, flatten::Any = 1.0, j_start::Int = 1, j_end::Int = size(dmatrix,2))

    costfun = cost_functional(cfspec)
    jmax = size(dmatrix,2)
    if ~(1 <= j_start <= j_end <= jmax)
        error("j_start and j_end should satisfy the inequality 1 <= j_start <= j_end <= "*string(jmax))
    end

    # Select the c2f best basis
    (dvecc2f, BSc2f) = ghwt_c2f_bestbasis(dmatrix, GP, cfspec = cfspec, flatten = flatten, j_start = j_start, j_end = j_end)
    costc2f = costfun(dvecc2f)

    # Select the f2c best basis
    (dvecf2c, BSf2c) = ghwt_f2c_bestbasis(dmatrix, GP, cfspec = cfspec, flatten = flatten, j_start = jmax + 1 - j_end, j_end = jmax + 1 - j_start)
    costf2c = costfun(dvecf2c)

    # Compare the coarse-to-fine and fine-to-coarse best-bases
    if costc2f >= costf2c
        dvec = dvecf2c
        BS = BSf2c
    else
        dvec = dvecc2f
        BS = BSc2f
    end

    # Return things
    return dvec, BS
end # of function ghwt_bestbasis

"""
    function GHWT_jkl(GP::GraphPart, drow::Int, dcol::Int; c2f::Bool = true)

Generate the (j,k,l) indices for the GHWT basis vector corresponding to the coefficient dmatrix(drow,dcol)

### Input Arguments
* `GP`: a GraphPart object
* `drow`: the row of the expansion coefficient
* `dcol`: the column of the expansion coefficient

### Output Argument
* `j`: the level index of the expansion coefficient
* `k`: the subregion index of the expansion coefficient
* `l`: the tag of the expansion coefficient
"""
function GHWT_jkl(GP::GraphPart, drow::Int, dcol::Int; c2f::Bool = true)
    if c2f
        j = dcol - 1

        k = findfirst(GP.rs[:,dcol] .> drow)
        k -= 2

        if isempty(GP.tag)
            ghwt_core!(GP)
        end
        l = GP.tag[drow,dcol]
    else
        jmax = size(GP.rs, 2)
        j = jmax - dcol

        IX = fine2coarse!(GP; indp = true)
        k = findfirst(GP.rs[:,j+1] .> IX[drow,dcol])
        k -= 2

        l = GP.tagf2c[drow,dcol]
    end
    return j,k,l
end # of function GHWT_jkl

end # of module GHWT
