"""
    natural_eigdist(๐ฝ, ๐, Q; ฮฑ = 1.0, T = :Inf,
                    input_format = :zero_measures, distance = :DAG,
                    edge_weight = 1, edge_length = 1)

compute natural distances between graph Laplacian eigenvectors.

# Input Arguments
- `๐ฝ::Matrix{Float64}`: matrix of (weighted) graph Laplacian eigenvectors.
- `๐::Vector{Float64}`: vector of eigenvalues.
- `Q::Matrix{Float64}`: unweighted incidence matrix of the graph.
- `ฮฑ::Float64`: ROT parameter. (default: `1.0`)
- `T::Any`: TSD parameter, i.e., the stopping time T in K_functional (default: `:Inf`)
- `input_format::Symbol`: options: `:zero_measures`, `:pmf1` and `:pmf2` (default: `:zero_measures`)
- `distance::Symbol`: options: `:ROT`, `:HAD`, `:DAG` and `:TSD` (default: `:DAG`)
- `edg_length::Any`: vector of edge lengths (default: 1 represents unweighted graphs)
- `edge_weight::Any`: the weights vector, which stores the affinity weight of
    each edge (default: `1` represents unweighted graphs).

# Output Argument
- `dis::Matrix{Float64}`: the distance matrix, dis[i,j] = d(๐แตขโโ, ๐โฑผโโ).

"""
function natural_eigdist(๐ฝ, ๐, Q; ฮฑ = 1.0, T = :Inf,
                         input_format = :zero_measures, distance = :DAG,
                         edge_weight = 1, edge_length = 1)
    N = size(Q, 1)
    P = deepcopy(๐ฝ)
    if input_format == :zero_measures
        P[:, 1] .= 0
    elseif input_format == :pmf1
        P = P.^2
    elseif input_format == :pmf2
        P = exp.(P) ./ sum(exp.(P), dims = 1)
    else
        @error("input_format does not support $(input_format)!")
        return
    end

    if distance == :ROT
        D = eigROT_Distance(P, Q; edge_length = edge_length, ฮฑ = ฮฑ)
    elseif distance == :HAD
        D = eigHAD_Distance(๐ฝ, ๐)
    elseif distance == :DAG
        D = eigDAG_Distance(๐ฝ, Q, N; edge_weight = edge_weight)
    elseif distance == :TSD
        tL = Q * Q'
        t๐, t๐ฝ = eigen(Matrix(tL))
        D = eigTSD_Distance(P, t๐ฝ, t๐, Q; length = edge_length, T = T)
    else
        error("distance does not support $(distance)!")
        return
    end

    return D
end
