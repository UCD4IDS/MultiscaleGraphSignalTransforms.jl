function rising_cutoff(t)
    r = 0
    if t <= -1
        r = 0
    elseif t >= 1
        r = 1
    else
        r = sin(π / 4 * (1 + sin(π / 2 * t)))
    end
    return r
end

function standardize_fiedler(v)
    v1 = deepcopy(v)
    v1 ./= norm(v1, 2)
    if v1[1] < 0
        v1 = -v1
    end
    return round.(v1; digits = 14)
end


function Lrw_eigenvec(W; nev = 6)
    # N = size(W, 1)
    # D = Diagonal(sum(W, dims = 1)[:])
    # try
    #     val, vtmp = eigs(D - W, D,
    #                      nev = nev, sigma = eps(), v0 = ones(N) / sqrt(N))
    #     vtmp = vtmp[:, sortperm(val)[2:end]]
    #     vtmp ./= sqrt.(sum(vtmp.^2; dims = 1))
    #     vtmp *= Diagonal(1 .- (vtmp[1, :] .< 0) .* 2)
    # catch
    #
    # return round.(vtmp; digits = 14)
    deg = sum(W, dims = 1)[:]  # weighted degree vector
    if minimum(deg) <= 10^3 * eps()
        L = Diagonal(deg) - W
        𝛌, 𝚽 = eigen(L)
        standardize_eigenvector_signs!(𝚽)
        return round.(𝚽[:, 2:nev]; digits = 14)
    end
    Lsym = Diagonal(deg.^(-1/2)) * (Diagonal(deg) - W) * Diagonal(deg.^(-1/2))
    𝛌sym, 𝚽sym = eigen(Lsym)
    𝚽rw = Diagonal(deg.^(-1/2)) * 𝚽sym
    𝚽rw ./= sqrt.(sum(𝚽rw.^2; dims = 1))
    𝚽rw *= Diagonal(1 .- (𝚽rw[1, :] .< 0) .* 2)
    return round.(𝚽rw[:, 2:nev]; digits = 14)
end

function sortnodes_inmargin(active_region, Na, v, W, idx; sign = :positive)
    ind = sortperm(v[active_region]; rev = (sign == :positive))
    # if there is a tie, use the more dims eigenmaps and sort by lexical order
    if length(unique(v[ind])) < length(active_region)
        emb = Lrw_eigenvec(W[idx, idx]; nev = 3)'
        # for the second dimension coordinate we sort
        # the pos region and the neg region in the same order
        emb_tp = [Tuple(emb[:, i] .* [1, (sign == :positive) * 2 - 1]) for i in 1:length(idx)]
        ind = sortperm(emb_tp[active_region]; rev = (sign == :positive))
    end
    return ind[(end - Na + 1):end]
end


function find_pairinds(W; ϵ::Float64 = 0.2, idx = 1:size(W, 1), used_node = Set())
    v = partition_fiedler(W[idx, idx])[2]
    v = standardize_fiedler(v)
    vmax = norm(v, Inf)
    N = length(v)
    pos_active_region = (Int64)[]
    neg_active_region = (Int64)[]

    tol = 10^3 * eps()
    tol *= ((sum(v .>= tol) > sum(v .<= -tol)) * 2 - 1)

    for i in 1:N
        if idx[i] ∈ used_node
            continue
        end
        if tol <= v[i] < ϵ * vmax
            push!(pos_active_region, i)
        end
        if -ϵ * vmax < v[i] < tol
            push!(neg_active_region, i)
        end
    end
    Np = length(pos_active_region)
    Nn = length(neg_active_region)
    Na = min(Nn, Np)  # number of pair inds in action region

    # indp = sortperm(v[pos_active_region]; rev = true)[(end - Na + 1):end]
    # indn = sortperm(v[neg_active_region])[(end - Na + 1):end]
    indp = sortnodes_inmargin(pos_active_region, Na, v, W, idx; sign = :positive)
    indn = sortnodes_inmargin(neg_active_region, Na, v, W, idx; sign = :negative)
    pair_inds = vcat(idx[pos_active_region[indp]]',
                     idx[neg_active_region[indn]]')

    return pair_inds, v
end

function pair_inds_shadding(W, GP; ϵ = 0.3, J = 1)
    rs = GP.rs
    inds = GP.inds
    (N, jmax) = Base.size(inds)

    used_node = Set()
    shading = zeros(N)
    for j = 1:J
        regioncount = count(!iszero, rs[:, j]) - 1
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            pair_inds = find_pairinds(W; ϵ = ϵ, idx = inds[indr, j], used_node = used_node)[1]
            for i in 1:size(pair_inds, 2)
                pv, nv = pair_inds[:, i]
                # if pv ∈ used_node || nv ∈ used_node
                #     continue
                # end
                if shading[pv] == 0
                    shading[pv] = J + 3 - j
                end
                if shading[nv] == 0
                    shading[nv] = -(J + 3 - j)
                end
            end
            union!(used_node, Set(pair_inds[:]))
        end
    end
    return shading
end



function keep_folding!(U, used_node, W, GP; ϵ = 0.3, j = 1)
    if j == 0
        return U
    end

    rs = GP.rs
    inds = GP.inds
    N = Base.size(inds, 1)

    regioncount = count(!iszero, rs[:, j]) - 1
    for r = 1:regioncount
        indr = rs[r, j]:(rs[r + 1, j] - 1)
        if length(indr) == 1
            continue
        end
        pair_inds, v = find_pairinds(W; ϵ = ϵ, idx = inds[indr, j], used_node = used_node)
        vmax = norm(v, Inf)
        for i in 1:size(pair_inds, 2)
            pv, nv = pair_inds[:, i]
            if pv ∈ used_node || nv ∈ used_node
                continue
            end
            # use the half distance between the 1D embeddings,
            # i.e., t = const⋅(v[pv] - v[nv] / 2), to compute the rising
            # cutoff function, which satisfy r(t)^2 + r(-t)^2 = 1
            t = (v[findfirst(inds[indr, j] .== pv)] - v[findfirst(inds[indr, j] .== nv)]) / (2 * ϵ * vmax)
            # t = v[findfirst(inds[indr, j] .== pv)] / (ϵ * vmax)
            U[pv, pv] = rising_cutoff(t)
            U[pv, nv] = rising_cutoff(-t)
            U[nv, pv] = -rising_cutoff(-t)
            U[nv, nv] = rising_cutoff(t)
        end
        union!(used_node, Set(pair_inds[:]))
    end
end


function unitary_folding_operator(W, GP; ϵ = 0.3, J = 1)
    rs = GP.rs
    inds = GP.inds
    (N, jmax) = Base.size(inds)

    U = Matrix{Float64}(I, N, N)
    used_node = Set()

    for j = 1:J
        keep_folding!(U, used_node, W, GP; ϵ = ϵ, j = j)
    end

    return U
end



"""
    lp_ngwp(𝚽::Matrix{Float64}, W_dual::SparseMatrixCSC{Float64, Int64},
            GP_dual::GraphPart; ϵ::Float64 = 0.3)

construct the lapped NGWP and GP.tag in place.

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors
- `W_dual::SparseMatrixCSC{Float64, Int64}`: weight matrix of the dual graph
- `GP_dual::GraphPart`: GraphPart object of the dual graph

# Output Argument
- `wavelet_packet::Array{Float64,3}`: the lapped NGWP dictionary. The first
    index is for selecting wavelets at a fixed level; the second index is for
    selecting the level `j`; the third index is for selecting elements in the
    wavelet vector.

"""
function lp_ngwp(𝚽::Matrix{Float64}, W_dual::SparseMatrixCSC{Float64, Int64},
                 GP_dual::GraphPart; ϵ::Float64 = 0.3)
    rs = GP_dual.rs
    inds = GP_dual.inds
    (N, jmax) = Base.size(inds)

    GP_dual.tag = zeros(Int, N, jmax)
    GP_dual.tag[:, 1] = Vector{Int}(0:(N - 1))

    Uf = Matrix{Float64}(I, N, N)
    used_node = Set()

    wavelet_packet = zeros(N, jmax, N)
    wavelet_packet[:, 1, :] = Matrix{Float64}(I, N, N)
    for j = 2:jmax
        regioncount = count(!iszero, rs[:, j]) - 1
        # Uf = unitary_folding_operator(W_dual, GP_dual; ϵ = ϵ, J = j - 1)
        keep_folding!(Uf, used_node, W_dual, GP_dual; ϵ = ϵ, j = j - 1)
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            GP_dual.tag[indr, j] = Vector{Int}(0:(length(indr) - 1))
            wavelet_packet[indr, j, :] = const_meyer_wavelets(𝚽, Uf; idx = inds[indr, j])'
        end
    end
    return wavelet_packet
end

function const_meyer_wavelets(𝚽, Uf; idx = 1:size(Uf, 1))
    N = size(𝚽, 1)
    # 1. nice interpretable way by smooth orthogonal projector
    # assemble smooth orthogonal projector custom to nodes `idx`
    P = Uf[idx, :]' * Uf[idx, idx]
    if diag(P) == χ(idx, N)
        B = 𝚽[:, idx]
    else
        # folding the eigenspace, i.e., 𝚽's column space
        Y = 𝚽 * P
        # find its column space's orthogonal basis
        B = svd(Y).U
    end
    # 2. alternative faster way by orthogonal unfolding operator
    # B = 𝚽 * Uf[:, idx]
    # perform varimax rotation to get the meyer_wavelets
    Wavelets = varimax(B)
    return Wavelets
end

"""
    lp_ngwp_analysis(G::GraphSig, 𝚽::Matrix{Float64}, W_dual::SparseMatrixCSC{Float64, Int64},
            GP_dual::GraphPart; ϵ::Float64 = 0.3)

perform the LP-NGWP transform of the graph signal(s) `G.f`.

# Input Arguments
- `G::GraphSig`: a GraphSig object
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors
- `W_dual::SparseMatrixCSC{Float64, Int64}`: weight matrix of the dual graph
- `GP_dual::GraphPart`: GraphPart object of the dual graph

# Output Argument
- `wavelet_packet::Array{Float64,3}`: the lapped NGWP dictionary. The first
    index is for selecting wavelets at a fixed level; the second index is for
    selecting the level `j`; the third index is for selecting elements in the
    wavelet vector.

"""
function lp_ngwp_analysis(G::GraphSig, 𝚽::Matrix{Float64}, W_dual::SparseMatrixCSC{Float64, Int64},
                 GP_dual::GraphPart; ϵ::Float64 = 0.3)
    rs = GP_dual.rs
    inds = GP_dual.inds
    (N, jmax) = Base.size(inds)
    fcols = size(G.f, 2)

    GP_dual.tag = zeros(Int, N, jmax)
    GP_dual.tag[:, 1] = Vector{Int}(0:(N - 1))

    Uf = Matrix{Float64}(I, N, N)
    used_node = Set()

    dmatrix = zeros(N, jmax, fcols)
    dmatrix[:, 1, :] = G.f
    for j = 2:jmax
        regioncount = count(!iszero, rs[:, j]) - 1
        keep_folding!(Uf, used_node, W_dual, GP_dual; ϵ = ϵ, j = j - 1)
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            GP_dual.tag[indr, j] = Vector{Int}(0:(length(indr) - 1))
            dmatrix[indr, j, :] = const_meyer_wavelets(𝚽, Uf; idx = inds[indr, j])' * G.f
        end
    end
    return dmatrix
end
