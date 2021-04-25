using NGWP, MTSG, Plots, LightGraphs, MultivariateStats
using LaTeXStrings

Nx, Ny = 7, 3
G = LightGraphs.grid([Nx, Ny]); N = nv(G);
L = Matrix(laplacian_matrix(G))
Q = incidence_matrix(G; oriented = true)
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽.*sign.(𝚽[1,:])'; # sign of DCT

∇𝚽 = Q' * 𝚽
W = 1.0 * adjacency_matrix(G)


grid2eig_ind = [1,2,3,6,8,12,15,4,5,7,9,13,16,18,10,11,14,17,19,20,21];
eig2grid_ind = sortperm(grid2eig_ind);
eig2dct = Array{Int64,3}(undef, Nx, Ny, 2);
for i = 1:Nx; for j = 1:Ny; eig2dct[i,j,1] = i-1; eig2dct[i,j,2] = j-1; end; end
eig2dct = reshape(eig2dct, (N, 2)); eig2dct = eig2dct[eig2grid_ind, :];


function grid7x3_mds_heatmaps(E, 𝚽; Nx = 7, Ny = 3, backend = :pyplot,
                              annotate_ind = 1:N, plotOrder = 1:N)
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

    # generate Grid7x3 2D MDS heatmaps plot
    if backend == :gr
        gr(dpi = 200)
    elseif backend == :pyplot
        pyplot(dpi = 200)
    elseif backend == :plotlyjs
        plotlyjs(dpi = 200)
    else
        @error("backend does not support $(backend)!")
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

# display(vcat((1:N)', eig2dct'))
