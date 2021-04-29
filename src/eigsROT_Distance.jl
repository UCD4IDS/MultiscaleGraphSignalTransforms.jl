"""
    eigsROT_Distance(P::Matrix{Float64}, W::SparseMatrixCSC{Float64, Int64}, X::Matrix{Float64}; α::Float64 = 1.0)

computes the sROT distance matrix of P's column vectors on a (unweighted) tree.

# Input Argument
- `P::Matrix{Float64}`: a matrix whose columns are vector measures with the same total mass.
- `W::SparseMatrixCSC{Float64, Int64}`: the weight matrix of the tree.
- `X::Matrix{Float64}`: the node positions (i-th row represents node `i`'s location)
- `α::Float64`: default is 1.0. ROT parameter.

# Output Argument
- `dist_sROT::Matrix{Float64}`: distance matrix, dist_sROT[i,j] = d_{sROT}(pᵢ, pⱼ; α).
- `Ws::SparseMatrixCSC{Float64, Int64}`: the weight matrix of the simplified tree.
- `Xs::Matrix{Float64}`: the node locations of the simplified tree
- `𝚯::Matrix{Float64}`: the shortened pmfs from `P`.

"""
function eigsROT_Distance(P::Matrix{Float64}, W::SparseMatrixCSC{Float64, Int64}, X::Matrix{Float64}; α::Float64 = 1.0)
    G = SimpleGraph(W)
    # check if the input weight matrix forms a tree
    if ne(G) + 1 != nv(G) || !is_connected(G)
        @error("input graph is not a tree.")
    end
    d = degree(G)
    # find index of junction nodes
    ijc = findall(d .> 2)
    # cut the graph into several disconnected subgraphs
    Wc = deepcopy(W)
    for i in ijc; Wc[i, :] .= 0; Wc[:, i] .= 0; end
    # find the indices for each subgraph
    Ind = find_subgraph_inds(Wc)
    # low dimensional pmfs
    𝚯 = Ind' * P
    # the centroids of subgraphs
    Xs = Diagonal(1 ./ sum(Ind, dims = 1)[:]) * Ind' * X
    # build Gs, i.e., the graph of subgraphs or the simplified tree
    Gs = Graph(size(Ind, 2))
    Ns = nv(Gs)
    # index of Gs's junction nodes
    ijcs = []
    for k = 1:Ns
        supportind = findall(Ind[:, k] .== 1)
        if issubset(supportind, ijc)
            push!(ijcs, k)
        end
    end
    # connect the branches to the junctions to form the simplified tree
    for k in setdiff(1:Ns, ijcs)
        supportind = findall(Ind[:,k] .== 1)
        for i = 1:length(ijc)
            if sum(W[supportind, ijc[i]]) > 0
                add_edge!(Gs, Edge(k, ijcs[i]))
            end
        end
    end
    Ws = 1.0 * adjacency_matrix(Gs)
    # compute the sROT distance matrix
    Qs = incidence_matrix(Gs; oriented = true)
    dist_sROT = eigROT_Distance(𝚯, Qs; edge_length = 1, α = α)
    return dist_sROT, Ws, Xs, 𝚯
end


"""
    find_subgraph_inds(Wc::SparseMatrixCSC{Float64, Int64})

find all subgraph indices of a tree. (subgraph includes: branches and junctions)

# Input Argument
- `Wc::SparseMatrixCSC{Float64, Int64}`: the weight matrix of the tree chopped by
    the junctions.

# Output Argument
- `Ind::Matrix{Float64}`: a matrix whose columns represent the subgraph node
    indices in binary form.

"""
function find_subgraph_inds(Wc::SparseMatrixCSC{Float64, Int64})
    N = size(Wc, 1)
    Dc = Diagonal(sum(Wc, dims=1)[:])
    Lc = Matrix(Dc - Wc)
    𝛌c, 𝚽c = eigen(Lc)
    p = findall(𝛌c .< 100 * eps())
    Uc = 𝚽c[:, p]
    # use heat diffusion at t = ∞ and the siginal `1:N` to figure out all branch indices
    a = round.(1e8 * Uc * Uc' * [k for k = 1:N])
    ind = unique(a)
    Ind = zeros(N, length(ind))
    for k = 1:length(ind); Ind[findall(a .== ind[k]), k] .= 1; end
    return Ind
end
