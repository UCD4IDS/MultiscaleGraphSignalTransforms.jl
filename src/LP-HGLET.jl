"""
    function LPHGLET_Synthesis(dvec::Vector{Float64}, GP::GraphPart, BS::BasisSpec, G::GraphSig; gltype::Symbol = :L, ϵ::Float64 = 0.3)

Perform Lapped-HGLET Synthesis transform

### Input Arguments
* `dvec`: the expansion coefficients corresponding to the chosen basis
* `GP`: a GraphPart object
* `BS`: a BasisSpec object
* `G`: a GraphSig object
* `gltype`: :L or :Lsym, indicating which eigenvectors are used
* `ϵ`: relative action bandwidth (default: 0.3)

### Output Argument
* `f`: the reconstructed signal
* `GS`: the reconstructed GraphSig object
"""
function LPHGLET_Synthesis(dvec::Matrix{Float64}, GP::GraphPart, BS::BasisSpec, G::GraphSig; gltype::Symbol = :L, ϵ::Float64 = 0.3)
    # Preliminaries
    W = G.W
    inds = GP.inds
    rs = GP.rs
    N = size(W, 1)
    jmax = size(rs, 2)
    Uf = Matrix{Float64}(I, N, N)
    used_node = Set()

    # fill in the appropriate entries of dmatrix
    dmatrix = dvec2dmatrix(dvec, GP, BS)

    f = zeros(size(dmatrix[:, jmax, :]))


    # Perform the synthesis transform
    for j = 1:jmax
        regioncount = count(!iszero, rs[:,j]) - 1
        # assemble orthogonal folding operator at level j - 1
        keep_folding!(Uf, used_node, W, GP; ϵ = ϵ, j = j - 1)
        for r = 1:regioncount
            # indices of current region
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            # indices of current region's nodes
            indrs = inds[indr, j]
            # number of nodes in current region
            n = length(indrs)

            # only proceed forward if coefficients do not exist
            if (j == jmax || count(!iszero, dmatrix[indr, j + 1, :]) == 0) && count(!iszero, dmatrix[indr, j, :]) > 0
                # compute the eigenvectors
                W_temp = W[indrs,indrs]
                D_temp = Diagonal(vec(sum(W_temp, dims = 1)))
                if gltype == :L
                    # compute the eigenvectors of L ==> svd(L)
                    v = svd(Matrix(D_temp - W_temp)).U
                elseif gltype == :Lsym
                    # check if one can assemble the Lsym
                    if minimum(sum(W[indrs, indrs], dims = 1)) > 10^3 * eps()
                        ### eigenvectors of L_sym ==> svd(L_sym)
                        D_temp_p = Diagonal(vec(sum(W_temp, dims = 1)).^(-0.5))
                        v = svd(Matrix(D_temp_p * (D_temp - W_temp) * D_temp_p)).U
                    else
                        ### eigenvectors of L ==> svd(L)
                        v = svd(Matrix(D_temp - W_temp)).U
                    end
                end
                v = v[:,end:-1:1]


                # standardize the eigenvector signs
                standardize_eigenvector_signs!(v)

                # construct unfolder operator custom to current region
                P = Uf[indrs, :]'

                # reconstruct the signal
                f += (P * v) * dmatrix[indr, j, :]

            end
        end
    end

    # creat a GraphSig object with the reconstructed data
    GS = deepcopy(G)
    replace_data!(GS, f)

    return f, GS
end



"""
    function LPHGLET_Analysis_All(G::GraphSig, GP::GraphPart; ϵ::Float64 = 0.3)

For a GraphSig object 'G', generate the 2 matrices of Lapped-HGLET expansion coefficients
corresponding to the eigenvectors of L and Lsym

### Input Arguments
* `G`:  a GraphSig object
* `GP`: a GraphPart object
* `ϵ`: relative action bandwidth (default: 0.3)

### Output Argument
* `dmatrixlH`:        the matrix of expansion coefficients for L
* `dmatrixlHsym`:     the matrix of expansion coefficients for Lsym
* `GP`:              a GraphPart object
"""
function LPHGLET_Analysis_All(G::GraphSig, GP::GraphPart; ϵ::Float64 = 0.3)
    # Preliminaries
    W = G.W
    inds = GP.inds
    rs = GP.rs
    N = size(W, 1)
    jmax = size(rs, 2)
    fcols = size(G.f, 2)
    Uf = Matrix{Float64}(I, N, N)
    used_node = Set()
    dmatrixlH = zeros(N, jmax, fcols)
    dmatrixlHsym = deepcopy(dmatrixlH)

    for j = 1:jmax
        regioncount = count(!iszero, rs[:,j]) - 1
        # assemble orthogonal folding operator at level j - 1
        keep_folding!(Uf, used_node, W, GP; ϵ = ϵ, j = j - 1)
        for r = 1:regioncount
            # indices of current region
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            # indices of current region's nodes
            indrs = inds[indr, j]
            # number of nodes in current region
            n = length(indrs)

            # compute the eigenvectors
            W_temp = W[indrs,indrs]
            D_temp = Diagonal(vec(sum(W_temp, dims = 1)))
            ## eigenvectors of L ==> svd(L)
            v = svd(Matrix(D_temp - W_temp)).U
            ## eigenvectors of L_sym ==> svd(L_sym)
            if minimum(sum(W[indrs, indrs], dims = 1)) > 10^3 * eps()
                ### eigenvectors of L_sym ==> svd(L_sym)
                D_temp_p = Diagonal(vec(sum(W_temp, dims = 1)).^(-0.5))
                v_sym = svd(Matrix(D_temp_p * (D_temp - W_temp) * D_temp_p)).U
            else
                ### eigenvectors of L ==> svd(L)
                v_sym = deepcopy(v)
            end

            # standardize the eigenvector signs
            v = v[:,end:-1:1]
            standardize_eigenvector_signs!(v)
            v_sym = v_sym[:,end:-1:1]
            standardize_eigenvector_signs!(v_sym)

            # construct unfolding operator custom to current region
            P = Uf[indrs, :]'
            # obtain the expansion coefficients
            dmatrixlH[indr, j, :] = (P * v)' * G.f
            dmatrixlHsym[indr, j, :] = (P * v_sym)' * G.f
        end
    end

    return dmatrixlH, dmatrixlHsym

end



function standardize_eigenvector_signs!(v)
    # standardize the eigenvector signs for HGLET (different with NGWPs)
    for col = 1:size(v, 2)
        row = 1
        standardized = false
        while !standardized
            if v[row, col] > 10^3 * eps()
                standardized = true
            elseif v[row,col] < -10^3 * eps()
                v[:, col] = -v[:, col]
            else
                row += 1
            end
        end
    end
end

"""
    HGLET_dictionary(GP::GraphPart, G::GraphSig; gltype::Symbol = :L)

assemble the whole HGLET dictionary

### Input Arguments
* `GP`: a GraphPart object
* `G`:  a GraphSig object
* `gltype`: `:L` or `:Lsym`

### Output Argument
* `dictionary`: the HGLET dictionary

"""
function HGLET_dictionary(GP::GraphPart, G::GraphSig; gltype::Symbol = :L)
    N = size(G.W, 1)
    jmax = size(GP.rs, 2)
    dictionary = zeros(N, jmax, N)
    for j = 1:jmax
        BS = BasisSpec(collect(enumerate(j * ones(Int, N))))
        dictionary[:, j, :] = HGLET_Synthesis(Matrix{Float64}(I, N, N), GP, BS, G; gltype = gltype)[1]'
    end
    return dictionary
end

"""
    LPHGLET_dictionary(GP::GraphPart, G::GraphSig; gltype::Symbol = :L, ϵ::Float64 = 0.3)

assemble the whole LP-HGLET dictionary

### Input Arguments
* `GP`: a GraphPart object
* `G`:  a GraphSig object
* `gltype`: `:L` or `:Lsym`
* `ϵ`: relative action bandwidth (default: 0.3)

### Output Argument
* `dictionary`: the LP-HGLET dictionary

"""
function LPHGLET_dictionary(GP::GraphPart, G::GraphSig; gltype::Symbol = :L, ϵ::Float64 = 0.3)
    N = size(G.W, 1)
    jmax = size(GP.rs, 2)
    dictionary = zeros(N, jmax, N)
    for j = 1:jmax
        BS = BasisSpec(collect(enumerate(j * ones(Int, N))))
        dictionary[:, j, :] = LPHGLET_Synthesis(Matrix{Float64}(I, N, N), GP, BS, G; gltype = gltype, ϵ = ϵ)[1]'
    end
    return dictionary
end
