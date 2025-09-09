"""
    eigROT_Distance(P, Q; edge_length = 1, α = 1.0)

computes the ROT distance matrix of P's column vectors on a graph.

# Input Argument
- `P::Matrix{Float64}`: a matrix whose columns are vector measures with the same total mass.
- `Q::SparseMatrixCSC{Int64,Int64}`: the unweighted incidence matrix of the graph.
- `edge_length::Any`: the length vector (default: 1 represents unweighted graphs)
- `α::Float64`: default is 1.0. ROT parameter.

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, dis[i,j] = d_{ROT}(pᵢ, pⱼ; α).

"""
function eigROT_Distance(P::Matrix{Float64}, Q::SparseMatrixCSC{Int64,Int64};
                            edge_length::Any = 1, α::Float64 = 1.0)
    n = size(P, 2)
    total_mass = sum(P, dims = 1)[:]
    if norm(total_mass - total_mass[1] * ones(n), Inf) > 10^4 * eps()
        @error("P's column measures do not share the same total mass.")
        return
    end
    dis = zeros(n,n)
    le2 = [edge_length; edge_length]
    Q2 = [Q -Q]
    m2 = size(Q2, 2)
    for i = 1:(n - 1), j = (i + 1):n
        f = P[:, i] - P[:, j]
        md = Model(optimizer_with_attributes(HiGHS.Optimizer, "LogLevel" => 0))
        @variable(md, w[1:m2] >= 0.0)
        edge_length == 1 ? @objective(md, Min, sum(w)) : @objective(md, Min, sum(w .* le2))
        @constraint(md, Q2 * w .== f)
        JuMP.optimize!(md)
        wt = abs.(JuMP.value.(w))
        dis[i,j] = edge_length == 1 ? norm(wt.^α, 1) : norm((wt.^α) .* le2, 1)
    end
    return dis + dis'
end


"""
    ROT_Distance(A, B, Q; edge_length = 1, α = 1.0)

computes the ROT distance matrix from A's column vectors to B's column vectors.
If A, B are vector inputs, then it also returns the cost value and the optimal
transport plan.

# Input Argument
- `A::Any`: a vector or matrix whose columns are initial probability measures.
- `B::Any`: a vector or matrix whose columns are terminal probability measures.
- `Q::SparseMatrixCSC{Int64,Int64}`: the unweighted incidence matrix of the graph.
- `edge_length::Any`: the length vector (default: `1` represents unweighted graphs)
- `α::Float64`: default is `1.0`. ROT parameter.
- `retPlan::Bool`: an indicator if return the optimal plan (default: `false`)

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, dis[i,j] = d_{ROT}(aᵢ, bⱼ; α).

"""
function ROT_Distance(A::Any, B::Any, Q::SparseMatrixCSC{Int64,Int64};
                        edge_length::Any = 1, α::Float64 = 1.0, retPlan::Bool = false)
    m = ndims(A) > 1 ? size(A, 2) : 1
    n = ndims(B) > 1 ? size(B, 2) : 1
    dis = zeros(m, n)
    le2 = [edge_length; edge_length]
    Q2 = [Q -Q]
    m2 = size(Q2, 2)
    for i = 1:m, j = 1:n
        f = (ndims(A) > 1 && ndims(B) > 1) ? B[:,j] - A[:,i] : B - A
        md = Model(optimizer_with_attributes(HiGHS.Optimizer, "LogLevel" => 0))
        @variable(md, w[1:m2] >= 0.0)
        edge_length == 1 ? @objective(md, Min, sum(w)) : @objective(md, Min, sum(w .* le2))
        @constraint(md, Q2 * w .== f)
        JuMP.optimize!(md)
        wt = abs.(JuMP.value.(w))
        dis[i, j] = edge_length == 1 ? norm(wt.^α, 1) : norm((wt.^α) .* le2, 1)
        if ndims(A) == 1 && ndims(B) == 1
            if retPlan
                return dis[1, 1], wt
            else
                return dis[1, 1]
            end
        end
    end
    return dis
end
