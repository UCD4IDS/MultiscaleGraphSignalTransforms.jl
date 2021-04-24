"""
    dualgraph(dist::Matrix{Float64}; method::Symbol = :inverse, σ::Float64 = 1.0)

build the dual graph's weight matrix based on the given non-trivial eigenvector
metric.

# Input Arguments
- `dist::Matrix{Float64}`: eigenvector distance matrix
- `method::Symbol`: default is by taking inverse of the distance between
    eigenvectors. Ways to build the dual graph edge weights. Option: `:inverse`,
    `:gaussian`.
- `σ::Float64`: default is `1.0`. Gaussian variance parameter.

# Output Argument
- `G_star::GraphSig`: A `GraphSig` object containing the weight matrix of the
    dual graph.

"""
function dualgraph(dist::Matrix{Float64}; method::Symbol = :inverse, σ::Float64 = 1.0)
    N = Base.size(dist, 1)
    W_star = zeros(N, N)
    if method == :inverse
        for i = 1:(N - 1), j = (i + 1):N
            W_star[i, j] = 1 / dist[i, j]
        end
    elseif method == :gaussian
        for i = 1:(N - 1), j = (i + 1):N
            W_star[i, j] = exp(-dist[i, j] / σ^2)
        end
    else
        error("method must be :inverse or :gaussian.")
    end
    W_star = W_star + W_star'
    G_star = GraphSig(sparse(W_star); name = "dual graph")
    return G_star
end
