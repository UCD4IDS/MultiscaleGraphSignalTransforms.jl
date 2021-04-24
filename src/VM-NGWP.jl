"""
    vm_ngwp(𝚽::Matrix{Float64}, GP_star::GraphPart)

construct varimax NGWP and GP_star.tag in place.

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `GP_star::GraphPart`: GraphPart object of the dual graph

# Output Argument
- `wavelet_packet::Array{Float64,3}`: the varimax NGWP. The first index is for
    selecting wavelets at a fixed level; the second index is for selecting the
    level `j`; the third index is for selecting elements in the wavelet vector.

"""
function vm_ngwp(𝚽::Matrix{Float64}, GP_star::GraphPart)
    rs = GP_star.rs
    inds = GP_star.inds
    (N, jmax) = Base.size(inds)

    GP_star.tag = zeros(Int, N, jmax)
    GP_star.tag[:, 1] = Vector{Int}(0:(N - 1))

    wavelet_packet = zeros(N, jmax, N)
    wavelet_packet[:, 1, :] = Matrix{Float64}(I, N, N)

    for j = 2:jmax
        regioncount = count(!iszero, rs[:, j]) - 1
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            GP_star.tag[indr, j] = Vector{Int}(0:(length(indr) - 1))
            wavelet_packet[indr, j, :] = varimax(𝚽[:, inds[indr, j]])'
        end
    end

    return wavelet_packet
end
