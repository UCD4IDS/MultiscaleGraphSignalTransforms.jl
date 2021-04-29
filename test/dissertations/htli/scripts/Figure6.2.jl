cd(@__DIR__); include("setups/grid7x3.jl");
pyplot(dpi = 200)

# DAG pseudo-metric
distDAG = eigDAG_Distance(𝚽, Q, N)

# MDS embedding into R²
D = distDAG
E = transform(fit(MDS, D, maxoutdim=2, distances=true))

# set up all heatmap plots' positions
dx = 0.01; dy = dx;
xej = zeros(Nx, N); yej=zeros(Ny, N);
a = 5.0; b = 9.0;
for k = 1:N
    xej[:,k] = LinRange(E[1,k] - Ny * a * dx, E[1, k] + Ny * a * dx, Nx)
    yej[:,k] = LinRange(E[2,k] - a * dy, E[2, k] + a * dy, Ny)
end

##
plot()
for k=1:N
    heatmap!(xej[:, k], yej[:, k], reshape(𝚽[:, k], (Nx, Ny))', c = :viridis,
                colorbar = false, ratio = 1, annotations = (xej[4, k],
                yej[3, k] + b*dy, text(latexstring("\\varphi_{",
                string(eig2dct[k, 1]), ",", string(eig2dct[k, 2]), "}"), 10)))
end
plt = plot!(xlim = [-1.4, 1.3], ylim = [-1.4, 1.3], grid = false, clims = (-0.4, 0.4))
# first level partition
p1x = [-0.2, 1.0, NaN]; p1y = [1.3, -1.0, NaN];
plot!(p1x, p1y, c = :red, legend = false, width = 3)
# second level partition
p2x = [-1.0, 0.2, NaN, 0.4, 1.2, NaN]; p2y = [-0.8, 0.45, NaN, 0.25, 0.2, NaN];
plot!(p2x, p2y, c=:orange, legend = false, width = 2)
plt = current()
savefig(plt, "../figs/Grid7x3_DAG_2levels_partition.png")
