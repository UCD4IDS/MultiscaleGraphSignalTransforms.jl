"""
    eigDAG_Distance(𝚽, Q, numEigs; edge_weights = 1)

compute DAG distances between pairwise graph Laplacian eigenvectors.

# Input Arguments
- `𝚽::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, 𝜙ⱼ₋₁ (j = 1,...,size(𝚽,1)).
- `Q::Matrix{Float64}`: incidence matrix of the graph.
- `numEigs::Int64`: number of eigenvectors considered.
- `edge_weight::Any`: default value is 1, stands for unweighted graph
    (i.e., all edge weights equal to 1). For weighted graph, edge_weight is the
    weights vector, which stores the affinity weight of each edge.

# Output Argument
- `dis::Matrix{Float64}`: a numEigs x numEigs distance matrix, dis[i,j] = d_DAG(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).

"""
function eigDAG_Distance(𝚽, Q, numEigs; edge_weight = 1)
    dis = zeros(numEigs, numEigs)
    abs_∇𝚽 = abs.(Q' * 𝚽)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_∇𝚽[:,i]-abs_∇𝚽[:,j],2) : norm((abs_∇𝚽[:,i]-abs_∇𝚽[:,j]).*sqrt.(edge_weight),2)
    end
    return dis + dis'
end

function eigDAG_Distance_normalized(𝚽,Q,numEigs; edge_weight = 1)
    dis = zeros(numEigs,numEigs)
    abs_∇𝚽 = abs.(Q' * 𝚽)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_∇𝚽[:,i]-abs_∇𝚽[:,j],2) : norm((abs_∇𝚽[:,i]-abs_∇𝚽[:,j]).*sqrt.(edge_weight),2)
        dis[i,j] /= norm(𝚽[:,i] .* 𝚽[:,j],2)
    end
    return dis + dis'
end
