using Graphs, SimpleWeightedGraphs, LinearAlgebra
"""
    SunFlowerGraph(; N = 400)

SUNFLOWERGRAPH construct a simple weighted sunflower graph with N vertices. Edge weights are the reciprocal of Euclidean distances.

# Input Arguments
- `N::Int64`: default is 400, the number of vertices. Requires N > 26.

# Output Argument
- `G::SimpleWeightedGraph{Int64,Float64}`: a simple weighted graph of the sunflower.
- `L::Matrix{Float64}`: the weighted unnormalized graph Laplacian matrix.
- `X::Matrix{Float64}`: a matrix whose i-th row represent the 2D coordinates of the i-th node.

"""
function SunFlowerGraph(; N = 400)
    c=1.0/N; θ=(sqrt(5.0)-1)*π;
    X = zeros(N,2); for k=1:N X[k,:]=c*(k-1)*[cos((k-1)*θ) sin((k-1)*θ)]; end

    G = SimpleWeightedGraph(N)
    for k = 2:8
        G.weights[1,k] = 1/norm(X[1,:] - X[k,:])
        G.weights[k,1] = 1/norm(X[1,:] - X[k,:])
    end

    for k = 1:N
        if k+8 <= N
            G.weights[k,k+8] = 1/norm(X[k,:] - X[k+8,:])
            G.weights[k+8,k] = 1/norm(X[k,:] - X[k+8,:])
        end
        if k+13 <= N
            G.weights[k,k+13] = 1/norm(X[k,:] - X[k+13,:])
            G.weights[k+13,k] = 1/norm(X[k,:] - X[k+13,:])
        end
    end
    W = weights(G) #weighted adjacency_matrix
    Lw = Diagonal(sum(W, dims = 2)[:]) - W #weighted laplacian_matrix
    return G, Lw, X
end
