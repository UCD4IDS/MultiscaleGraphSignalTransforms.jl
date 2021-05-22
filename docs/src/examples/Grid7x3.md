# Metrics of Graph Laplacian Eigenvectors on ``P_7 \times P_3``

## Set up
```@example grid
using MultiscaleGraphSignalTransforms, LightGraphs, MultivariateStats
using Plots, LaTeXStrings, LinearAlgebra; pyplot(dpi = 200)

# compute the graph related quantities
Nx, Ny = 7, 3
G = LightGraphs.grid([Nx, Ny]); N = nv(G);
L = Matrix(laplacian_matrix(G))
Q = incidence_matrix(G; oriented = true)
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽.*sign.(𝚽[1,:])'; # sign of DCT
∇𝚽 = Q' * 𝚽
W = 1.0 * adjacency_matrix(G)

# manually set up the mapping between 1D ordering and 2D ordering
grid2eig_ind = [1,2,3,6,8,12,15,4,5,7,9,13,16,18,10,11,14,17,19,20,21];
eig2grid_ind = sortperm(grid2eig_ind);
eig2dct = Array{Int64,3}(undef, Nx, Ny, 2);
for i = 1:Nx; for j = 1:Ny; eig2dct[i,j,1] = i-1; eig2dct[i,j,2] = j-1; end; end
eig2dct = reshape(eig2dct, (N, 2)); eig2dct = eig2dct[eig2grid_ind, :];

# show the comparison between 1D ordering vs. 2D ordering of the eigenvectors
## 1D ordering
plot(layout = Plots.grid(3, 7), title = "Non-decreasing eigenvalue ordering")
for i in 1:N
    heatmap!(reshape(𝚽[:, i], (Nx, Ny))', c = :viridis, cbar = false,
                clims = (-0.4,0.4), frame = :none, ratio = 1, ylim = [0, Ny + 1],
                title = latexstring("\\phi_{", i-1, "}"), titlefont = 12,
                subplot = i)
end
plot!(size = (849, 350)) # hide
savefig("grid_1d_order.svg"); nothing # hide

## 2D ordering
plot(layout = Plots.grid(3, 7), title = "Natural frequency ordering")
for i in 1:N
    k = grid2eig_ind[i]
    heatmap!(reshape(𝚽[:,k], (Nx, Ny))', c = :viridis, cbar = false,
                clims = (-0.4,0.4), frame = :none, ratio = 1, ylim = [0, Ny + 1],
                title = latexstring("\\varphi_{", string(eig2dct[k,1]),
                ",", string(eig2dct[k,2]), "}"), titlefont = 12, subplot = i)
end
plot!(size = (849, 350)) # hide
savefig("grid_2d_order.svg"); nothing # hide
```
![](grid_1d_order.svg)

![](grid_2d_order.svg)




```@example grid
function grid7x3_mds_heatmaps(E, 𝚽; Nx = 7, Ny = 3, annotate_ind = 1:N, plotOrder = 1:N)
    # set up all heatmap plots' positions
    max_x = maximum(E[1, :]); min_x = minimum(E[1, :])
    width_x = max_x - min_x
    max_y = maximum(E[2, :]); min_y = minimum(E[2, :])
    width_y = max_y - min_y
    dx = 0.005 * width_x; dy = dx;
    xej = zeros(Nx, N); yej=zeros(Ny, N);
    a = 5.0; b = 7.0;
    for k = 1:N
        xej[:,k] = LinRange(E[1,k] - Ny * a * dx, E[1, k] + Ny * a * dx, Nx)
        yej[:,k] = LinRange(E[2,k] - a * dy, E[2, k] + a * dy, Ny)
    end

    plot()
    for k in plotOrder
        if k in annotate_ind
            heatmap!(xej[:, k], yej[:, k], reshape(𝚽[:, k], (Nx, Ny))', c = :viridis,
                     colorbar = false, ratio = 1, annotations = (xej[4, k],
                     yej[3, k] + b*dy, text(latexstring("\\varphi_{",
                     string(eig2dct[k, 1]), ",", string(eig2dct[k, 2]), "}"), 10)))
        else
            heatmap!(xej[:, k], yej[:, k], reshape(𝚽[:, k], (Nx, Ny))', c = :viridis,
                     colorbar = false, ratio = 1)
        end
    end
    plt = plot!(xlim = [min_x - 0.12 * width_x, max_x + 0.12 * width_x],
                ylim = [min_y - 0.16 * width_y, max_y + 0.16 * width_y],
                grid = false, clims = (-0.4, 0.4),
                xlab = "X₁", ylab = "X₂")
    return plt
end
nothing # hide
```

## ROT distance
Before we measure the ROT distance between the eigenvectors, we convert them to probability mass functions by taking entrywise squares.
After we got the ROT distance matrix of the eigenvectors, we visualize the arrangement of the eigenvectors in ``\mathbb{R}^{2}`` via [the classical MDS embedding](https://en.wikipedia.org/wiki/Multidimensional_scaling#Classical_multidimensional_scaling).
```@example grid
## ROT distance
D = natural_eigdist(𝚽, 𝛌, Q; α = 0.5, input_format = :pmf1, distance = :ROT)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
grid7x3_mds_heatmaps(E, 𝚽)
plot!(size = (815, 611)) # hide
```

## DAG distance
We organize the eigenvectors by the DAG distance.
```@example grid
D = natural_eigdist(𝚽, 𝛌, Q; distance = :DAG)
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
grid7x3_mds_heatmaps(E, 𝚽)
plot!(size = (815, 611)) # hide
```

## TSD distance
We organize the eigenvectors by the TSD distance with the parameter ``T = 0.1``.
```@example grid
D = natural_eigdist(𝚽, 𝛌, Q; T = 0.1, distance = :TSD)  # T = 0.1
E = transform(fit(MDS, D, maxoutdim=2, distances=true))
grid7x3_mds_heatmaps(E, 𝚽)
plot!(size = (815, 611)) # hide
```
