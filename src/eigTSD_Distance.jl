"""
    eigTSD_Distance(P::Matrix{Float64}, 𝚽::Matrix{Float64}, 𝛌::Vector{Float64},
                    Q::SparseMatrixCSC{Int64,Int64}; length::Any = 1,
                    T::Any = :Inf, tol::Float64 = 1e-5)

computes the TSD distance matrix of P's column vectors on a graph.

# Input Argument
- `P::Matrix{Float64}`: vector measures with the same total mass 0.
- `𝚽::Matrix{Float64}`: matrix of the unweighted graph Laplacian eigenvectors.
- `𝛌::Vector{Float64}`: vector of eigenvalues.
- `Q::SparseMatrixCSC{Int64,Int64}`: the unweighted incidence matrix.
- `length::Any`: vector of edge lengths (default: `1` represents unweighted graphs)
- `T::Any`: the stopping time T in K_functional (default: `:Inf`)
- `tol::Float64`: tolerance for integral convergence (default: `1e-5`)

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, d_{TSD}(φᵢ, φⱼ; T).

"""
function eigTSD_Distance(P::Matrix{Float64}, 𝚽::Matrix{Float64}, 𝛌::Vector{Float64},
                         Q::SparseMatrixCSC{Int64,Int64}; length::Any = 1,
                         T::Any = :Inf, tol::Float64 = 1e-5)
    N, ncols = Base.size(P)
    total_mass = sum(P, dims = 1)[:]
    if norm(total_mass - total_mass[1] * ones(ncols), Inf) > 10^4 * eps()
        @error("P's column measures do not share the same total mass.")
        return
    end
    # initialize the distance matrix
    dis = zeros(ncols, ncols)
    # store gradient of 𝚽 to avoid repeated computation
    ∇𝚽 = Q' * 𝚽

    for i = 1:(ncols - 1), j = (i + 1):ncols
        dis[i, j] = K_functional(P[:, i], P[:, j], 𝚽, ∇𝚽, 𝛌; length = length,
                                    T = T, tol = tol)[1]
    end
    return dis + dis'
end

"""
    K_functional(𝐩::Vector{Float64}, 𝐪::Vector{Float64}, 𝚽::Matrix{Float64},
                 ∇𝚽::Matrix{Float64}, 𝛌::Vector{Float64}; length::Any = 1,
                 T::Any = :Inf, dt::Float64 = 0.5/maximum(𝛌), tol::Float64 = 1e-5)

computes the K_functional between two vector meassures 𝐩 and 𝐪 on a graph.

# Input Argument
- `𝐩::Vector{Float64}`: the source vector measure.
- `𝐪::Vector{Float64}`: the destination vector measure.
- `𝚽::Matrix{Float64}`: matrix of the unweighted graph Laplacian eigenvectors.
- `∇𝚽::Matrix{Float64}`: gradient of unweighted graph Laplacian eigenvectors.
- `𝛌::Vector{Float64}`: vector of eigenvalues.
- `length::Any`: vector of edge lengths (default: 1 represents unweighted graphs)
- `T::Any`: the stopping time T in K_functional (default: :Inf)
- `tol::Float64`: tolerance for convergence (default: 1e-5)

# Output Argument
- `K::Float64`: TSD distance d_{TSD}(p, q; T).
- `E::Float64`: an estimated upper bound on the absolute error. In general,
    `E <= tol * norm(K)`.

"""
function K_functional(𝐩::Vector{Float64}, 𝐪::Vector{Float64}, 𝚽::Matrix{Float64},
                        ∇𝚽::Matrix{Float64}, 𝛌::Vector{Float64}; length::Any = 1,
                        T::Any = :Inf, tol::Float64 = 1e-5)
    if abs(sum(𝐩 - 𝐪)) > 10^4 * eps()
        @error("𝐩 and 𝐪 do not have the same total mass.")
    end
    f₀ = 𝐪 - 𝐩
    # store b to avoid repeated computation
    b = 𝚽' * f₀

    if T == :Inf
        # choose stopping time T such that exp(-λ₁T) = tol
        T = -log(tol) / 𝛌[2]
    end
    # use QuadGK.jl to numerically evaluate the integral
    K, E = quadgk(t -> weighted_1norm(∇f(t, ∇𝚽, b, 𝛌), length), 0, T, rtol=tol)
    return K, E
end

function ∇f(t, ∇𝚽, b, 𝛌)
    gu = ∇𝚽 * (exp.(-t * 𝛌) .* b)
    return gu
end

function weighted_1norm(x, length)
    if length == 1
        return norm(x, 1)
    else
        return sum(abs.(x) .* length)
    end
end
