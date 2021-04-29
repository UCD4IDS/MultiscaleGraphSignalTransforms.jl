using MultiscaleGraphSignalTransforms, LightGraphs, Plots, LaTeXStrings

## Build Graph
N = 128; G = path_graph(N)
X = zeros(N,2); X[:, 1] = 1:N
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L); standardize_eigenvectors!(𝚽)
W = 1.0 * adjacency_matrix(G)

G_Sig = GraphSig(W, xy = X)
GP = partition_tree_fiedler(G_Sig; swapRegion = false)

function anti_diag(A)
    N = size(A, 1)
    return [A[i, N+1-i] for i = 1:N]
end


## construct 1D smooth orthogonal projector
ϵa = 16
pair_inds = vcat(vec((64 + ϵa):-1:65)', vec((64 - ϵa + 1):64)')
Uf = Matrix{Float64}(I, N, N)
β = 64.5

for i in 1:size(pair_inds, 2)
    pv, nv = pair_inds[:, i]
    t = abs(pv - β) / ϵa
    Uf[pv, pv] = rising_cutoff(t)
    Uf[pv, nv] = rising_cutoff(-t)
    Uf[nv, pv] = -rising_cutoff(-t)
    Uf[nv, nv] = rising_cutoff(t)
end

P0 = Uf' * diagm(χ(1:64, N)) * Uf
P1 = Uf' * diagm(χ(65:N, N)) * Uf


##
distDCT = zeros(N,N)
for i in 1:N-1, j = i+1:N
    distDCT[i,j] = abs(i-j)
end
distDCT = distDCT + distDCT'


function path_spectrogram(f, D, 𝚽; c = 0.01)
    N = length(f)
    dmatrix = zeros(N, N)
    for l = 1:N
        P = 𝚽 * diagm(nat_spec_filter(l, D; σ = c * maximum(D))) * 𝚽'
        for x = 1:N
            ψ = P * spike(x, N)
            ψ ./= norm(ψ, 2)
            dmatrix[l, x] = ψ' * f
        end
    end
    heatmap(abs.(dmatrix), c = :thermal, ratio = 1, xlim = [0.5, N+1],
            ylim = [0.5, N], xlab = latexstring("x"), ylab = latexstring("l"),
            guidefontsize = 18)
    xticks!([16:16:128;], [string(k) for k in 16:16:128])
    yticks!([1;16:16:128], vcat("DC", [string(k) for k in 15:16:127]))
    plt = plot!(size = (500, 400))
    return plt
end
