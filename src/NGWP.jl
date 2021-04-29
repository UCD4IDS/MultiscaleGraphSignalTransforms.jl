
"""
    ngwp_analysis(G::GraphSig, wavelet_packet::Array{Float64,3})

For a GraphSig object `G`, generate the matrix of NGWP expansion coefficients.

# Input Arguments
- `G::GraphSig`: an input GraphSig object
- `wavelet_packet::Array{Float64,3}`: the varimax wavelets packet.

# Output Argument
- `dmatrix::Array{Float64,3}`: the expansion coefficients matrix.

"""
function ngwp_analysis(G::GraphSig, wavelet_packet::Array{Float64,3})
    f = G.f
    fcols = size(f, 2)
    (N, jmax, ) = Base.size(wavelet_packet)

    dmatrix = zeros(N, jmax, fcols)
    dmatrix[:, 1, :] = f

    for j = 2:jmax
        for i = 1:N
            dmatrix[i, j, :] = f' * wavelet_packet[i, j, :]
        end
    end

    return dmatrix
end


"""
    const_proj_wavelets(𝚽,vlist,elist; method = "Modified Gram-Schmidt with Lp Pivoting")

construct projection wavelets, i.e., project δ ∈ vlist onto span({φⱼ| j ∈ elist}).

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `vlist::Array{Int}`: the list of considered node indices.
- `elist::Array{Int}`: the list of considered eigenvector indices.
- `method::Symbol`: default is `:MGSLp`. other options: `:IP` (Iterative-Projection),
    `GS` (Gram Schmidt).

# Output Argument
- `Wav::Matrix{Float64}`: a matrix whose columns are projected wavelet vectors.

"""
function const_proj_wavelets(𝚽, vlist, elist; method = :MGSLp)
    if length(vlist) == 1
        return 𝚽[:, elist]
    end
    N = size(𝚽, 1)
    m = length(vlist)
    Wav = zeros(N, m)

    B = 𝚽[:, elist]

    if method == :IP
        for k in 1:length(vlist)
            wavelet = Proj(spike(vlist[k], N), B)
            Wav[:, k] .= wavelet ./ norm(wavelet)
            B = wavelet_perp_Matrix(wavelet, B)
        end
    elseif method == :GS
        P = B * B'
        for k in 1:length(vlist)
            wavelet = P * spike(vlist[k], N)
            Wav[:,k] .= wavelet ./ norm(wavelet)
        end
        Wav, complement_dim = gram_schmidt(Wav)
        if complement_dim != 0
            complement_space = B * nullspace(Wav' * B)
            Wav = hcat(Wav, complement_space)
        end
    elseif method == :MGSLp
        P = B * B'
        for k in 1:length(vlist)
            wavelet = P * spike(vlist[k], N)
            Wav[:,k] .= wavelet ./ norm(wavelet)
        end
        Wav, complement_dim = mgslp(Wav)
        if complement_dim != 0
            complement_space = B * nullspace(Wav' * B)
            Wav = hcat(Wav, complement_space)
        end
    else
        error("Do not support method ", method)
    end

    return Wav
end

"""
    function NGWP_jkl(GP_star::GraphPart, drow::Int, dcol::Int)

Generate the (j,k,l) indices for the NGWP basis vector corresponding to the coefficient dmatrix(drow,dcol)

### Input Arguments
* `GP_star::GraphPart`: a GraphPart object of the dual grpah
* `drow::Int`: the row of the expansion coefficient
* `dcol::Int`: the column of the expansion coefficient

### Output Argument
* `j`: the level index of the expansion coefficient
* `k`: the subregion in dual graph's index of the expansion coefficient
* `l`: the tag of the expansion coefficient
"""
function NGWP_jkl(GP_star::GraphPart, drow::Int, dcol::Int)
    (j, k, l) = GHWT_jkl(GP_star, drow, dcol)
    return j,k,l
end

"""
    (dvec_ngwp, BS_ngwp) = ngwp_bestbasis(dmatrix::Array{Float64,3}, GP_star::GraphPart;
                                          cfspec::Any = 1.0, flatten::Any = 1.0,
                                          j_start::Int = 1, j_end::Int = size(dmatrix, 2),
                                          useParent::Bool = true)

Select the best basis from the matrix of NGWP expansion coefficients.

# Input Arguments
- `dmatrix::Array{Float64,3}`: the matrix of expansion coefficients
- `GP_star::GraphPart`: an input GraphPart object of the dual graph
- `cfspec::Any`: the specification of cost functional to be used (default = 1.0,
    i.e., 1-norm)
- `flatten::Any`: the method for flattening vector-valued data to scalar-valued
    data (default = 1.0, i.e, 1-norm)
- `useParent::Bool`: the flag to indicate if we update the selected best basis
    subspace to the parent when parent and child have the same cost (default = false)

#  Output Arguments
- `dvec_ngwp::Matrix{Float64}`: the vector of expansion coefficients corresponding
    to the NGWP best basis
- `BS_ngwp::BasisSpec`: a BasisSpec object which specifies the NGWP best basis
"""
function ngwp_bestbasis(dmatrix::Array{Float64,3}, GP_star::GraphPart;
                            cfspec::Any = 1.0, flatten::Any = 1.0,
                            j_start::Int = 1, j_end::Int = size(dmatrix, 2),
                            useParent::Bool = false)
    dvec_ngwp, BS_ngwp = ghwt_c2f_bestbasis(dmatrix, GP_star; cfspec = cfspec,
                            flatten = flatten, j_start = j_start, j_end = j_end,
                            useParent = useParent)
    BS_ngwp.description = "NGWP Best Basis"
    return dvec_ngwp, BS_ngwp
end
