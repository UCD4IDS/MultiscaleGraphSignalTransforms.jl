"""
    findminimum(v, n)

FINDMINIMUM finds the first n smallest elements' indices.

# Input Arguments
- `v::Array{Float64}`: the candidate values for selection.
- `n::Int`: number of smallest elements for consideration.

# Output Argument
- `idx::Array{Int}`: n smallest elements' indices.

"""
function findminimum(v, n)
    idx = sortperm(v)[1:n]
    return idx
end


"""
    spike(i, N)

SPIKE gives the N-dim spike vector with i-th element equals 1.

# Input Arguments
- `i::Int`: index for one.
- `N::Int`: dimension of the target spike vector.

# Output Argument
- `v::Array{Float64}`: the N-dim spike vector with i-th element equals 1.

"""
function spike(i, N)
    v = zeros(N)
    v[i] = 1
    return v
end

"""
    characteristic(list, N)

CHARACTERISTIC gives the characteristic function in n-dim vector space with
    values of index in list equal to 1.

# Input Arguments
- `list::Array{Int}`: list of indices.
- `N::Int`: dimension of the target vector.

# Output Argument
- `v::Array{Float64}`: the n-dim characteristic vector with values of index in list equal to 1.

"""
function characteristic(list, N)
    v = zeros(N)
    v[list] .= 1.0
    return v
end

function χ(list, N)
    v = zeros(N)
    v[list] .= 1.0
    return v
end


"""
    heat_sol(f0,𝚽,Σ,t)

HEAT\\_SOL gives the solution of heat partial differential equation with initial condition u(⋅, 0) = f0

# Input Arguments
- `f0::Array{Float64}`: initial condition vector.
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors, served as graph Fourier transform matrix
- `Σ::Array{Int}`: diagonal matrix of eigenvalues.
- `t::Float`: time elapse.

# Output Argument
- `u::Array{Float64}`: the solution vector at time t

"""
function heat_sol(f0, 𝚽, Σ, t)
    u = 𝚽 * (exp.(-t .* Σ) .* 𝚽' * f0)
    return u
end


"""
    freq_band_matrix(ls,n)

FREQ\\_BAND\\_MATRIX provides characteristic diagonal matrix, which is useful for spectral graph filters design.

# Input Arguments
- `ls::Array{Int}`: list of indices.
- `n::Int`: dimension of the target vector.

# Output Argument
- `D::Array{Float64}`: the zero/one diagonal matrix.

"""
function freq_band_matrix(ls, n)
    f = characteristic(list, n)
    return Diagonal(f)
end


"""
    scatter_gplot(X; marker = nothing, ms = 4, plotOrder = :normal, c = :viridis, subplot = 1)

SCATTER\\_GPLOT generates a scatter plot figure, which is for quick viewing of a graph signal.
SCATTER\\_GPLOT!(X; ...) adds a plot to `current` one.

# Input Arguments
- `X::Matrix{Float64}`: points locations, can be 2-dim or 3-dim.
- `marker::Array{Float64}`: default is `nothing`. Present different colors given
    different signal value at each node.
- `ms::Array{Float64}`: default is `4`. Present different node sizes given
    different signal value at each node.
- `shape::Symbol`: default is `:none`. Shape of the markers.
- `mswidth::Number`: default is `0`. Width of the marker stroke border.
- `msalpha::Number`: default is `nothing`. The opacity for the marker stroke.
- `plotOrder::Symbol`: default is `:normal`. Optional choices `:s2l` or `:l2s`, i.e.,
    plots from the smallest value of `marker` to the largest value or the other way around.
- `c::Symbol`: default is `:viridis`. Colors of the markers.
- `subgplot::Int`: default is `1`. The subplot index.

"""
function scatter_gplot(X; marker = nothing, ms = 4, shape = :none, mswidth = 0,
                          msalpha = nothing, plotOrder = :normal, c = :viridis, subplot = 1)
    N, dim = size(X)
    if marker != nothing
        if size(marker) == (N,) || size(marker) == (N, 1)
            marker = marker[:]  # reshape N x 1 matrix to a vector of length N
        else
            error("marker only accepts a vector of length $(N) or a matrix of size $(N) x 1.")
        end
        if plotOrder == :normal
            idx = 1:N
        elseif plotOrder == :s2l
            idx = sortperm(marker)
        elseif plotOrder == :l2s
            idx = sortperm(marker, rev = true)
        else
            error("plotOrder only supports for :normal, :s2l, or :l2s.")
        end
        X = X[idx, :]
        marker = marker[idx]
        if length(ms) > 1
            ms = ms[idx]
        end
    end
    if dim == 2
        scatter(X[:, 1], X[:, 2], marker_z = marker, ms = ms, shape = shape, c = c,
                legend = false, mswidth = mswidth, msalpha = msalpha, cbar = true,
                aspect_ratio = 1, grid = false, subplot = subplot)
    elseif dim == 3
        scatter(X[:, 1], X[:, 2], X[:, 3], marker_z = marker, ms = ms, shape = shape, c = c,
                legend = false, mswidth = mswidth, msalpha = msalpha, cbar = true,
                aspect_ratio = 1, grid = false, subplot = subplot)
    else
        error("Dimension Error: scatter_gplot only supports for 2-dim or 3-dim scatter plots.")
    end
end

function scatter_gplot!(X; marker = nothing, ms = 4, shape = :none, mswidth = 0,
                          msalpha = nothing, plotOrder = :normal, c = :viridis, subplot = 1)
    N, dim = size(X)
    if marker != nothing
        if size(marker) == (N,) || size(marker) == (N, 1)
            marker = marker[:]  # reshape N x 1 matrix to a vector of length N
        else
            error("marker only accepts a vector of length $(N) or a matrix of size $(N) x 1.")
        end
        if plotOrder == :normal
            idx = 1:N
        elseif plotOrder == :s2l
            idx = sortperm(marker)
        elseif plotOrder == :l2s
            idx = sortperm(marker, rev = true)
        else
            error("plotOrder only supports for :normal, :s2l, or :l2s.")
        end
        X = X[idx, :]
        marker = marker[idx]
        if length(ms) > 1
            ms = ms[idx]
        end
    end
    if dim == 2
        scatter!(X[:, 1], X[:, 2], marker_z = marker, ms = ms, shape = shape, c = c,
                legend = false, mswidth = mswidth, msalpha = msalpha, cbar = true,
                aspect_ratio = 1, grid = false, subplot = subplot)
    elseif dim == 3
        scatter!(X[:, 1], X[:, 2], X[:, 3], marker_z = marker, ms = ms, shape = shape, c = c,
                legend = false, mswidth = mswidth, msalpha = msalpha, cbar = true,
                aspect_ratio = 1, grid = false, subplot = subplot)
    else
        error("Dimension Error: scatter_gplot! only supports for 2-dim or 3-dim scatter plots.")
    end
end



"""
    cat_plot(X; marker = nothing, ms = 4)

CAT\\_PLOT generates a scatter plot figure for cat example, which is for quick
viewing of a graph signal within a specific range (i.e., xlims, ylims, zlims).
CAT\\_PLOT!(X; ...) adds a plot to `current` one.

# Input Arguments
- `X::Matrix{Float64}`: 3-dim points.
- `marker::Array{Float64}`: default is nothing. Present different colors given
    different signal value at each node.
- `ms::Array{Float64}`: default is 4. Present different node sizes given
    different signal value at each node.

"""
function cat_plot(X; marker = nothing, ms = 4)
    scatter(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = :viridis,
        legend = false, cbar = true, aspect_ratio = 1, xlims = [-100, 100],
        ylims = [-100, 100], zlims = [-100, 100])
end

function cat_plot!(X; marker = nothing, ms = 4)
    scatter!(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = :viridis,
        legend = false, cbar = true, aspect_ratio = 1, xlims = [-100, 100],
        ylims = [-100, 100], zlims = [-100, 100])
end


"""
    approx_error_plot(DVEC::Vector{Vector{Float64}}; frac::Float64 = 0.50)

draw relative approx. error w.r.t. fraction of coefficients retained (FCR)

# Input Arguments
- `DVEC::Vector{Vector{Float64}}`: a list of expansion coefficients.
- `fraction::Float64`: default is 0.5.

"""
function approx_error_plot(DVEC::Vector{Vector{Float64}}; frac::Float64 = 0.50)
    plot(xaxis = "Fraction of Coefficients Retained",
            yaxis = "Relative Approximation Error")
    T = ["Laplacian", "HGLET", "Haar", "Walsh", "GHWT_c2f", "GHWT_f2c", "eGHWT",
            "PC-NGWP", "VM-NGWP"]
    L = [(:dashdot, :red), (:solid, :brown), (:dashdot, :orange),
            (:dashdot, :pink), (:solid, :gray), (:solid,:green),
            (:solid,:blue), (:solid,:purple), (:solid,:black)]
    for i = 1:length(DVEC)
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = max.(sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm, 1e-12)
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac),
                label = T[i], line = L[i], linewidth = 2, grid = false)
    end
end


"""
    getall_expansioncoeffs(G_Sig::GraphSig, GP_star::GraphPart, VM_NGWP::Array{Float64,3}, PC_NGWP::Array{Float64,3}, 𝚽::Matrix{Float64})

get all expansion coefficients of `f` via all methods in NGWP.jl and MTSG.jl

# Input Arguments
- `G_Sig::GraphSig`: GraphSig of the primal graph
- `GP_star::GraphPart`: GraphPart of the dual graph
- `VM_NGWP::Array{Float64,3}`: varimax NGWP.
- `PC_NGWP::Array{Float64,3}`: pair-clustering NGWP.
- `𝚽::Matrix{Float64}`: the matrix of graph Laplacian eigenvectors.

"""
function getall_expansioncoeffs(G_Sig::GraphSig, GP_star::GraphPart, VM_NGWP::Array{Float64,3}, PC_NGWP::Array{Float64,3}, 𝚽::Matrix{Float64})
    ############# VM_NGWP
    dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
    dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_star)
    ############# PC_NGWP
    dmatrix_PC = ngwp_analysis(G_Sig, PC_NGWP)
    dvec_pc_ngwp, BS_pc_ngwp = ngwp_bestbasis(dmatrix_PC, GP_star)
    ############# Laplacian
    dvec_Laplacian = 𝚽' * G_Sig.f
    ############# plain Laplacian eigenvectors coefficients
    GP = partition_tree_fiedler(G_Sig)
    dmatrixH = HGLET_Analysis_All(G_Sig, GP)[1]
    dvec_hglet, BS_hglet, _ = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixH, costfun = 1)
    dmatrix = ghwt_analysis!(G_Sig, GP = GP)
    ############# Haar
    BS_haar = bs_haar(GP)
    dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)
    ############# Walsh
    BS_walsh = bs_walsh(GP)
    dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)
    ############# GHWT_c2f
    dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)
    ############# GHWT_f2c
    dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)
    ############# eGHWT
    dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)
    DVEC = [dvec_Laplacian[:], dvec_hglet[:], dvec_haar[:], dvec_walsh[:],
            dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:], dvec_pc_ngwp[:],
            dvec_vm_ngwp[:]]
    return DVEC
end

using Clustering
"""
    spectral_clustering(𝚽, M)

SPECTRAL_CLUSTERING return M graph clusters, i.e., {Vₖ| k = 1,2,...,M}.

# Input Argument
- `𝚽::Matrix{Float64}`: the matrix of graph Laplacian eigenvectors.
- `M::Int`: the number of graph clusters.

# Output Argument
- `clusters::Vector{Vector{Int}}`: graph cluster indices.

"""
function spectral_clustering(𝚽, M)
    if M < 2
        return [1:size(𝚽, 1)]
    end
    cluster_indices = assignments(kmeans(𝚽[:, 2:M]', M))
    clusters = Vector{Vector{Int}}[]
    for k in 1:M
        push!(clusters, findall(cluster_indices .== k)[:])
    end
    return clusters
end

"""
    transform2D(X; s = 1, t = [0,0])

TRANSFORM2D dilate each point of `X` by scale s and translate by 2D vector t.
"""
function transform2D(X; s = 1, t = [0,0])
    X1 = X .* s
    X2 = zeros(size(X))
    for i in 1:size(X,1)
        X2[i,1] = X1[i,1] + t[1]
        X2[i,2] = X1[i,2] + t[2]
    end
    return X2
end

"""
    NN_rendering(X, Img_Mat)

NN\\_RENDERING generates a rendering signal at each point of `X` from the image
`Img_Mat` by nearest neighbor method.
"""
function NN_rendering(X, Img_Mat)
    N = size(X,1)
    f = zeros(N)
    for i in 1:N
        nn_x, nn_y = Int(round(X[i, 2])), Int(round(X[i, 1]))
        if nn_x < 1 || nn_x > size(Img_Mat, 2) || nn_y < 1 || nn_y > size(Img_Mat, 1)
            print("Error: pixel out of boundary!")
            return
        end
        f[i] = Img_Mat[nn_x, nn_y]
    end
    return f
end

"""
    Bilinear_rendering(X, Img_Mat)

NN\\_RENDERING generates a rendering signal at each point of `X` from the image
`Img_Mat` by bilinear interpolation method.
"""
function Bilinear_rendering(X, Img_Mat)
    N = size(X,1)
    f = zeros(N)
    for i in 1:N
        x1, x2, y1, y2 = Int(floor(X[i, 2])), Int(floor(X[i, 2])) + 1, Int(floor(X[i, 1])), Int(floor(X[i, 1])) + 1
        x, y = X[i,2], X[i,1]
        F = [Img_Mat[x1,y1] Img_Mat[x1,y2]
            Img_Mat[x2,y1] Img_Mat[x2,y2]]
        prod_res = 1/((x2 - x1) * (y2 - y1)) * [x2-x x-x1] * F * [y2-y y-y1]'
        f[i] = prod_res[1,1]
    end
    return f
end

"""
    dct1d(k, N)

DCT1D returns k-th 1D DCT basis vector in Rᴺ.

# Input Arguments
- `k::Int64`: ord of DCT basis vector. k = 1,2,...,N.
- `N::Int64`: vector dimension.

# Output Argument
- `φ::Array{Float64}`: k-th 1D DCT basis vector in Rᴺ. (k is 1-indexed)
"""
function dct1d(k, N)
    φ = [cos(π * (k - 1) * (l + 0.5) / N) for l = 0:(N - 1)]
    return φ ./ norm(φ, 2)
end

"""
    dct2d_basis(N1, N2)

DCT2D\\_BASIS returns 2D DCT basis vectors in [0,1] x [0,1] with N1-1 and N2-1 subintervals respectively.

# Input Arguments
- `N1::Int64`: number of nodes in x-axis.
- `N2::Int64`: number of nodes in y-axis.

# Output Argument
- `𝚽::Matrix{Float64}`: 2D DCT basis vectors.
"""
function dct2d_basis(N1, N2)
    N = N1 * N2
    𝚽 = zeros(N, N)
    ind = 1
    for i in 1:N1, j in 1:N2
        φ₁, φ₂ = dct1d(i, N1), dct1d(j, N2)
        φ = reshape(φ₁ * φ₂', N)
        𝚽[:, ind] = φ
        ind += 1
    end
    return 𝚽
end

"""
    alternating_numbers(n)

ALTERNATING\\_NUMBERS e.g., n = 5, returns [1,5,2,4,3]; n = 6, returns [1,6,2,5,3,4]

# Input Arguments
- `n::Int64`: number of nodes in x-axis.

# Output Argument
- `arr::Array{Int64}`: result array.
"""
function alternating_numbers(n)
    mid = Int(ceil(n / 2))
    arr1 = 1:mid
    arr2 = n:-1:(mid + 1)
    arr = Array{Int64}(zeros(n))
    p1, p2 = 1, 1
    for i = 1:n
        if i % 2 == 1
            arr[i] = arr1[p1]
            p1 += 1
        else
            arr[i] = arr2[p2]
            p2 += 1
        end
    end
    return arr
end

"""
    compute_SNR(f, g)

COMPUTE\\_SNR, g = f + ϵ, SNR = 20 * log10(norm(f)/norm(g-f)).

# Input Arguments
- `f::Array{Float64}`: original signal.
- `g::Array{Float64}`: noisy signal.

# Output Argument
- `SNR::Float64`: SNR value.
"""
function compute_SNR(f, g)
    SNR = 20 * log10(norm(f) / norm(g - f))
    return SNR
end

"""
    sort_wavelets(A; onlyByLoc = false)

sort A's column wavelet vectors based on their focused nodes' indices. flip
signs via the cross correlation.

# Input Arguments
- `A::Matrix{Float64}`: whose column vectors are wavelets.

# Output Argument
- `A::Matrix{Float64}`: a matrix with sorted and sign flipped column.
"""
function sort_wavelets(A; onlyByLoc = false)
    ord = findmax(abs.(A), dims = 1)[2][:]
    idx = sortperm([j[1] for j in ord])
    A = A[:, idx]
    if onlyByLoc
        return A
    end

    N = size(A, 1)
    sgn = ones(size(A, 2))
    mid = Int(round(size(A, 2) / 2))
    A[:, mid] .*= (maximum(A[:, mid]) > -minimum(A[:, mid])) * 2 - 1
    for i in 1:size(A, 2)
        cor_res = crosscor(A[:, mid], A[:, i], -Int(ceil(N / 2)):Int(floor(N / 2)))
        if maximum(cor_res) < -minimum(cor_res)
            sgn[i] = -1
        end
    end
    A = A * Diagonal(sgn)
    return A
end

"""
    standardize_eigenvectors!(𝚽::Matrix{Float64})

standardize the signs of the eigenvectors such that 1) the term with the largest
magnitude is positive; 2) if the max == -min, then make sure the first non-zero
entry of the eigenvector is positive.

# Input Arguments
- `𝚽::Matrix{Float64}`: matrix of graph Laplacian eigenvectors

"""
function standardize_eigenvectors!(𝚽::Matrix{Float64})
    N, nev = size(𝚽)
    tol = 10^3 * eps()
    for l in 1:nev
        𝛟max = maximum(𝚽[:, l])
        𝛟min = minimum(𝚽[:, l])
        if abs(𝛟max + 𝛟min) < tol
            row = 1
            standardized = false
            while !standardized
                if 𝚽[row, l] > tol
                    standardized = true
                elseif 𝚽[row, l] < -tol
                    𝚽[:, l] = -𝚽[:, l]
                else
                    row += 1
                end
            end
        elseif 𝛟max < -𝛟min
            𝚽[:, l] *= -1
        end
    end
end


"""
getall_expansioncoeffs2(G_Sig::GraphSig, GP_star::GraphPart,
                                 GP_star_Lsym::GraphPart,
                                 VM_NGWP::Array{Float64,3},
                                 PC_NGWP::Array{Float64,3},
                                 LP_NGWP::Array{Float64,3},
                                 VM_NGWP_Lsym::Array{Float64,3},
                                 PC_NGWP_Lsym::Array{Float64,3},
                                 LP_NGWP_Lsym::Array{Float64,3},
                                 𝚽::Matrix{Float64},
                                 𝚽sym::Matrix{Float64})

get all expansion coefficients of `f` in dissertation via all methods in NGWP.jl
and MTSG.jl

# Input Arguments
- `G_Sig::GraphSig`: GraphSig of the primal graph
- `GP_star::GraphPart`: GraphPart of the dual graph
- `VM_NGWP::Array{Float64,3}`: varimax NGWP.
- `PC_NGWP::Array{Float64,3}`: pair-clustering NGWP.
- `LP_NGWP::Array{Float64,3}`: lapped NGWP.
- `VM_NGWP_Lsym::Array{Float64,3}`: Lsym version of varimax NGWP.
- `PC_NGWP_Lsym::Array{Float64,3}`: Lsym version of pair-clustering NGWP.
- `LP_NGWP_Lsym::Array{Float64,3}`: Lsym version of lapped NGWP.
- `𝚽::Matrix{Float64}`: the matrix of `L` eigenvectors.
- `𝚽sym::Matrix{Float64}`: the matrix of `Lsym` eigenvectors.

"""
function getall_expansioncoeffs2(G_Sig::GraphSig, GP_star::GraphPart,
                                 GP_star_Lsym::GraphPart,
                                 VM_NGWP::Array{Float64,3},
                                 PC_NGWP::Array{Float64,3},
                                 LP_NGWP::Array{Float64,3},
                                 VM_NGWP_Lsym::Array{Float64,3},
                                 PC_NGWP_Lsym::Array{Float64,3},
                                 LP_NGWP_Lsym::Array{Float64,3},
                                 𝚽::Matrix{Float64},
                                 𝚽sym::Matrix{Float64})
    ############# VM_NGWP
    dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
    dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_star)
    ############# PC_NGWP
    dmatrix_PC = ngwp_analysis(G_Sig, PC_NGWP)
    dvec_pc_ngwp, BS_pc_ngwp = ngwp_bestbasis(dmatrix_PC, GP_star)
    ############# LP_NGWP
    dmatrix_LP = ngwp_analysis(G_Sig, LP_NGWP)
    dvec_lp_ngwp, BS_lp_ngwp = ngwp_bestbasis(dmatrix_LP, GP_star)
    ############# VM_NGWP_Lsym
    dmatrix_VM_Lsym = ngwp_analysis(G_Sig, VM_NGWP_Lsym)
    dvec_vm_ngwp_Lsym, BS_vm_ngwp_Lsym = ngwp_bestbasis(dmatrix_VM_Lsym, GP_star_Lsym)
    ############# PC_NGWP_Lsym
    dmatrix_PC_Lsym = ngwp_analysis(G_Sig, PC_NGWP_Lsym)
    dvec_pc_ngwp_Lsym, BS_pc_ngwp_Lsym = ngwp_bestbasis(dmatrix_PC_Lsym, GP_star_Lsym)
    ############# LP_NGWP
    dmatrix_LP_Lsym = ngwp_analysis(G_Sig, LP_NGWP_Lsym)
    dvec_lp_ngwp_Lsym, BS_lp_ngwp_Lsym = ngwp_bestbasis(dmatrix_LP_Lsym, GP_star_Lsym)
    ############# Laplacian L
    dvec_Laplacian = 𝚽' * G_Sig.f
    ############# Laplacian Lsym
    dvec_Laplacian_sym = 𝚽sym' * G_Sig.f
    ############# HGLET
    GP = partition_tree_fiedler(G_Sig)
    dmatrixH, _, dmatrixHsym = HGLET_Analysis_All(G_Sig, GP)
    dvec_hglet, BS_hglet, trans_hglet = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixH, dmatrixHsym = dmatrixHsym, costfun = 1)
    ############# LP-HGLET
    dmatrixsH, dmatrixsHsym = LPHGLET_Analysis_All(G_Sig, GP; ϵ = 0.3)
    dvec_lphglet, BS_lphglet, trans_lphglet = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixsH, dmatrixHsym = dmatrixsHsym, costfun = 1)
    ############# GHWT dictionaries
    dmatrix = ghwt_analysis!(G_Sig, GP = GP)
    ############# Haar
    BS_haar = bs_haar(GP)
    dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)
    ############# Walsh
    BS_walsh = bs_walsh(GP)
    dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)
    ############# GHWT_c2f
    dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)
    ############# GHWT_f2c
    dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)
    ############# eGHWT
    dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)
    DVEC = [dvec_Laplacian[:], dvec_Laplacian_sym[:], dvec_hglet[:], dvec_lphglet[:],
            dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:],
            dvec_pc_ngwp[:], dvec_pc_ngwp_Lsym[:],
            dvec_vm_ngwp[:], dvec_vm_ngwp_Lsym[:],
            dvec_lp_ngwp[:], dvec_lp_ngwp_Lsym[:]]
    return DVEC
end

"""
    approx_error_plot2(DVEC::Vector{Vector{Float64}}; frac::Float64 = 0.50)

draw relative approx. error w.r.t. fraction of coefficients retained (FCR)
in dissertation

# Input Arguments
- `DVEC::Vector{Vector{Float64}}`: a list of expansion coefficients.
- `fraction::Float64`: default is 0.5.

"""
function approx_error_plot2(DVEC::Vector{Vector{Float64}}; frac::Float64 = 0.50)
    plot(xaxis = "Fraction of Coefficients Retained",
            yaxis = "Relative Approximation Error", size = (600, 500))
    T = ["eigenbasis-L", "eigenbasis-Lsym", "HGLET", "LP-HGLET", "Haar", "Walsh",
         "GHWT_c2f", "GHWT_f2c", "eGHWT", "PC-NGWP", "PC-NGWP-Lsym",
         "VM-NGWP", "VM-NGWP-Lsym", "LP-NGWP", "LP-NGWP-Lsym"]
    L = [(:dot, :red), (:dot, :magenta), (:solid, :teal), (:dashdot, :teal),
         (:dashdot, :orange), (:dashdot, :pink), (:solid, :gray),
         (:solid, :green), (:solid, :blue), (:solid, :purple), (:dash, :purple),
         (:solid, :black), (:dash, :black), (:solid, :orange), (:dash, :orange)]
    for i = 1:length(DVEC)
        if i in [5, 6, 7, 8]
            continue
        end
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = max.(sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm, 1e-12)
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac),
                label = T[i], line = L[i], linewidth = 2, grid = false)
    end
end
