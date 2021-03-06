"""
    eigDAG_Distance(๐ฝ, Q, numEigs; edge_weights = 1)

compute DAG distances between pairwise graph Laplacian eigenvectors.

# Input Arguments
- `๐ฝ::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, ๐โฑผโโ (j = 1,...,size(๐ฝ,1)).
- `Q::Matrix{Float64}`: incidence matrix of the graph.
- `numEigs::Int64`: number of eigenvectors considered.
- `edge_weight::Any`: default value is 1, stands for unweighted graph
    (i.e., all edge weights equal to 1). For weighted graph, edge_weight is the
    weights vector, which stores the affinity weight of each edge.

# Output Argument
- `dis::Matrix{Float64}`: a numEigs x numEigs distance matrix, dis[i,j] = d_DAG(๐แตขโโ, ๐โฑผโโ).

"""
function eigDAG_Distance(๐ฝ, Q, numEigs; edge_weight = 1)
    dis = zeros(numEigs, numEigs)
    abs_โ๐ฝ = abs.(Q' * ๐ฝ)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_โ๐ฝ[:,i]-abs_โ๐ฝ[:,j],2) : norm((abs_โ๐ฝ[:,i]-abs_โ๐ฝ[:,j]).*sqrt.(edge_weight),2)
    end
    return dis + dis'
end

function eigDAG_Distance_normalized(๐ฝ,Q,numEigs; edge_weight = 1)
    dis = zeros(numEigs,numEigs)
    abs_โ๐ฝ = abs.(Q' * ๐ฝ)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_โ๐ฝ[:,i]-abs_โ๐ฝ[:,j],2) : norm((abs_โ๐ฝ[:,i]-abs_โ๐ฝ[:,j]).*sqrt.(edge_weight),2)
        dis[i,j] /= norm(๐ฝ[:,i] .* ๐ฝ[:,j],2)
    end
    return dis + dis'
end
