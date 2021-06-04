"""
    nat_spec_filter(l, D; σ = 0.25 * maximum(D), method = :regular, thres = 0.2)

assemble the natural spectral graph filter centered at the l-th eigenvector via
the distance matrix `D`.

# Input Arguments
- `l::Int64`: index of the centered eigenvector
- `D::Matrix{Float64}`: non-trivial distance matrix of the eigenvectors
- `σ::Float64`: Gaussian window width parameter (default: `0.25 * maximum(D)`)
- `method::Symbol`: `:regular` or `:reduced` (default: `:regular`)
- `thres::Float64`: cutoff threshold ∈ (0, 1).

# Output Argument
- `𝛍::Vector{Float64}`: the natural spectral graph filter

"""
function nat_spec_filter(l, D; σ = 0.25 * maximum(D), method = :regular, thres = 0.2)
    d = D[:, l]
    𝛍 = exp.(-(d ./ σ).^2)
    𝛍 ./= sum(𝛍)
    if method == :reduced
        𝛍 .*= (𝛍 .> (thres * maximum(𝛍)))
    end
    return 𝛍
end

"""
    ngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D))

assemble the whole NGWF dictionary.

# Input Arguments
- `D::Matrix{Float64}`: non-trivial distance matrix of the eigenvectors
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors
- `σ::Float64`: Gaussian window width parameter (default: `0.25 * maximum(D)`)

# Output Argument
- `𝓤::Matrix{Float64}`: the NGWF dictionary

"""
function ngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D))
    N = size(D, 1)
    𝓤 = zeros(N, 0)
    for l = 1:N
        𝛍 = nat_spec_filter(l, D; σ = σ)
        P = 𝚽 * Diagonal(𝛍) * 𝚽'
        𝓤 = hcat(𝓤, P)
    end
    return 𝓤
end

"""
    rngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D), thres = 0.2)

assemble the reduced NGWF (rNGWF) dictionary.

# Input Arguments
- `D::Matrix{Float64}`: non-trivial distance matrix of the eigenvectors
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors
- `σ::Float64`: Gaussian window width parameter (default: `0.25 * maximum(D)`)
- `thres::Float64`: cutoff threshold ∈ (0, 1).

# Output Argument
- `𝓤::Matrix{Float64}`: the rNGWF dictionary
- `dic_l2x::Dict`: a dictionary to store the filtered locations by QR at the l-th
    centered eigenvector

"""
function rngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D), thres = 0.2)
    N = size(D, 1)
    𝓤 = zeros(N, 0)
    dic_l2x = Dict()
    for l = 1:N
        𝛍 = nat_spec_filter(l, D; σ = σ, method = :reduced, thres = thres)
        P = 𝚽 * Diagonal(𝛍) * 𝚽'
        dic_l2x[l] = qr(P, Val(true)).p[1:rank(P, rtol = 10^4 * eps())]
        𝓤 = hcat(𝓤, P[:, dic_l2x[l]])
    end
    return 𝓤, dic_l2x
end

function ngwf_vector(D, l, x, 𝚽; σ = 0.1 * maximum(D))
    N = size(𝚽, 1)
    P = 𝚽 * Diagonal(nat_spec_filter(l, D; σ = σ)) * 𝚽'
    ψ = P * spike(x, N)
    ψ ./= norm(ψ, 2)
    return ψ
end

"""
    frame_approx(f, U, V; num_kept = length(f))

approximate signal `f` by the frame `U`.

# Input Arguments
- `f::Vector{Float64}`: input graph signal
- `U::Matrix{Float64}`: a frame operator (matrix or dictionary)
- `V::Matrix{Float64}`: the dual frame operator of `U`
- `num_kept::Int64`: number of kept coefficients (NCR)

# Output Argument
- `rel_error::Vector{Float64}`: the relative errors
- `f_approx::Vector{Float64}`: the approximated signal

"""
function frame_approx(f, U, V; num_kept = length(f))
    g = U' * f
    ind = sortperm(g.^2; rev = true)[1:num_kept]
    f_approx = zeros(length(f))
    rel_error = [1.0]
    for γ in ind
        f_approx += g[γ] * V[:, γ]
        f_res = f - f_approx
        push!(rel_error, norm(f_res)/norm(f))
    end
    return rel_error, f_approx
end

"""
    rngwf_lx(dic_l2x)

find the sequential subindices of rNGWF vectors.

# Input Arguments
- `dic_l2x::Dict`: a dictionary to store the filtered locations by QR at the l-th
    centered eigenvector

# Output Argument
- `Γ::Vector{Tuple{Int64,Int64}}`: the sequential subindices of rNGWF vectors.

"""
function rngwf_lx(dic_l2x)
    N = length(dic_l2x)
    Γ = (Tuple{Int64,Int64})[]
    for l = 1:N
        for x in dic_l2x[l]
            push!(Γ, (l - 1, x))
        end
    end
    return Γ
end
