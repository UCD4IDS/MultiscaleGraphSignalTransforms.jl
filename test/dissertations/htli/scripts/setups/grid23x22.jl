using MLDatasets, LightGraphs, Plots, LaTeXStrings, MultiscaleGraphSignalTransforms
# Load local module
push!(LOAD_PATH, @__DIR__)
using pSGWT

example, label = MNIST.traindata(Float64, 28)
digit_img = example[4:26, 26:-1:5]
# heatmap(digit_img', ratio = 1, c=:viridis, frame = :none, xlim = [1, 22])

Nx, Ny = size(digit_img)
G = LightGraphs.grid([Nx, Ny]); N = nv(G);
W = Matrix(adjacency_matrix(G))
L = Matrix(laplacian_matrix(G))
Q = incidence_matrix(G; oriented = true)
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1, :])';

D = natural_eigdist(𝚽, 𝛌, Q; distance = :DAG)

function grid_eigenvector_plot(l)
    heatmap(reshape(𝚽[:, l], (Nx, Ny))', c = :viridis, ratio = 1, frame = :none,
        xlim = [1, Nx], size = (500, 400), cbar = true)
end

function grid_NGWFvector_plot(l, x, D, 𝚽; σ = 0.2 * maximum(D))
    heatmap(reshape(ngwf_vector(D, l, x, 𝚽; σ = σ)', (Nx, Ny))',
        c = :viridis, ratio = 1, frame = :none, xlim = [1, Nx], size = (500, 400))
end

function plot_edge!(A, B; style = :solid, subplot = 1, c = :red)
    plot!([A[1], B[1], NaN], [A[2], B[2], NaN], c = c, legend = false,
        width = 10, style = style, subplot = subplot)
end

function plot_square!(Nx, Ny; subplot = 1, c = :red)
    plot_edge!([-0.33, 0], [Nx + 1.24, 0]; subplot = subplot, c = c)
    plot_edge!([0, 0], [0, Ny + 1]; subplot = subplot, c = c)
    plot_edge!([-0.33, Ny + 1], [Nx + 1.24, Ny + 1]; subplot = subplot, c = c)
    plot_edge!([Nx + 0.92, 0], [Nx + 0.92, Ny + 1]; subplot = subplot, c = c)
end

function grid_vector_plot!(l, i, VECs)
    v = deepcopy(VECs[:, l])
    v ./= norm(v, 1)
    heatmap!(reshape(v, (Nx, Ny))', c = :viridis, ratio = 1, frame = :none,
        xlim = [-0.5, Nx+1.5], cbar = false, subplot = i)
end

function grid_vec_heatmap(VEC, Nx, Ny; l = 1)
    v = VEC[:, l]
    heatmap(reshape(v, (Nx, Ny))', c = :viridis, ratio = 1, frame = :none,
        xlim = [1, Nx], size = (500, 400), cbar = true)
end
