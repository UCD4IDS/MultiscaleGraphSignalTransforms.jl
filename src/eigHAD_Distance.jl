"""
    eigHAD_Distance(๐ฝ, ๐; indexEigs = 1:size(๐ฝ,2))

compute HAD "distance" (not really a distance) between pairwise graph Laplacian
eigenvectors, i.e., d_HAD(๐แตขโโ, ๐โฑผโโ) = โ(1 - a_HAD(๐แตขโโ, ๐โฑผโโ)ยฒ).

# Input Arguments
- `๐ฝ::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, ๐โฑผโโ (j = 1,...,size(๐ฝ,1)).
- `๐::Array{Float64}`: array of eigenvalues. (ascending order)
- `indexEigs::Int`: default is all eigenvectors, indices of eigenvectors considered.

# Output Argument
- `dis::Matrix{Float64}`: the HAD distance matrix, dis[i,j] = d_HAD(๐แตขโโ, ๐โฑผโโ).
"""
function eigHAD_Distance(๐ฝ, ๐; indexEigs = 1:size(๐ฝ,2))
    n = length(indexEigs)
    A = eigHAD_Affinity(๐ฝ, ๐; indexEigs = indexEigs)
    dis = sqrt.(ones(n, n) - A.^2)
    dis[diagind(dis)] .= 0
    return dis
end

function eigHAD_Distance_neglog(๐ฝ, ๐; indexEigs = 1:size(๐ฝ,2))
    A = eigHAD_Affinity(๐ฝ, ๐; indexEigs = indexEigs)
    n = size(A,1)
    dis = zeros(n,n)
    for i = 1:n, j = 1:n
        if A[i,j] == 0
            dis[i,j] = 1e9
        else
            dis[i,j] = -log(A[i,j])
        end
    end
    return dis
end

"""
    eigHAD_Affinity(๐ฝ, ๐; indexEigs = 1:size(๐ฝ,2))

EIGHAD_AFFINITY compute Hadamard (HAD) affinity between pairwise graph Laplacian eigenvectors.

# Input Arguments
- `๐ฝ::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, ๐โฑผโโ (j = 1,...,size(๐ฝ,1)).
- `๐::Array{Float64}`: array of eigenvalues. (ascending order)
- `indexEigs::Int`: default is all eigenvectors, indices of eigenvectors considered.

# Output Argument
- `A::Matrix{Float64}`: a numEigs x numEigs affinity matrix, A[i,j] = a_HAD(๐แตขโโ, ๐โฑผโโ).
"""
function eigHAD_Affinity(๐ฝ, ๐; indexEigs = 1:size(๐ฝ,2))
    N, numEigs = size(๐ฝ,1), length(indexEigs)
    indNoDC = setdiff(indexEigs, 1) # get rid of DC component
    J = length(indNoDC)
    A = zeros(J, J)
    for a in 1:J, b in a:J
        i, j = indNoDC[a], indNoDC[b]
        hadamardProd = ๐ฝ[:,i] .* ๐ฝ[:,j]
        if norm(hadamardProd,2) < 0.01/sqrt(N)
            continue
        end
        ฮป, ฮผ = ๐[i], ๐[j]
        xโ = 1 ./ (max(ฮป, ฮผ))
        # Find minimizer t
        result = optimize(t -> abs(exp(-t[1]*ฮป) + exp(-t[1]*ฮผ) - 1), [xโ], BFGS());
        t = Optim.minimizer(result)[1]
        # Compute Hadamard affinity
        heatEvolution = ๐ฝ * Diagonal(exp.(-t .* ๐)) * ๐ฝ' * hadamardProd
        A[a,b] = norm(heatEvolution,2) / (norm(hadamardProd,2) + 1e-6)
    end
    A = A + A'; for i in 1:J; A[i,i] /= 2; end

    if 1 in indexEigs
        # Set affinity measure of ๐โ with itself to be the maximum and equals to 1.
        A = matrix_augmentation(A)
        A[1,1] = maximum(A)
    end
    return A ./ maximum(A)
end

function matrix_augmentation(A)
    m, n = size(A)
    B = zeros(m+1, n+1)
    B[2:end, 2:end] = A
    return B
end
