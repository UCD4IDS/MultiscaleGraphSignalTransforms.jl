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



## Build Dual Graph
Gstar_Sig = dualgraph(distDAG)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual)

VM_NGWP = vm_ngwp(𝚽, GP_dual)

## level 2 VM-NGWP vectors
j = 3; W_VM = VM_NGWP[:, j, :]'

wav_kl = [[0 0];[0 1];[0 2];[1 0];[1 1];[1 2];[1 3];[2 0];[2 1];[2 2];[3 0];
            [3 1];[3 2];[3 3];[2 3];[2 4];[2 5];[2 6];[3 4];[3 5];[3 6]];
wav_kl = wav_kl[eig2grid_ind,:];

# reorder_ind = [2,3,1,5,7,4,6, 16,17,15,  9,11,8,10, 18,20,21,19, 13,14,12]
reorder_ind = [1,3,2,5,7,4,6, 8,10,9,17,15,18,16, 11,13,14,12,20,21,19]
W_VM = W_VM[:,reorder_ind[eig2grid_ind]];
sgn = ones(N); sgn[grid2eig_ind[[4,6,8,11,12,14]]] .= -1; W_VM = W_VM * Diagonal(sgn);

#################### Fig. 6
plot()
for k=1:N
    heatmap!(xej[:, k], yej[:, k], reshape(W_VM[:,k], (Nx, Ny))', c = :viridis,
        colorbar = false, ratio = 1, annotations = (xej[4, k], yej[3, k] + b * dy,
        text(latexstring("\\psi_{", string(wav_kl[k,1]), ",", string(wav_kl[k,2]), "}"), 10)))
end
plot!(aspect_ratio = 1, xlim = [-1.4, 1.3], ylim = [-1.4, 1.3], grid = false, clims=(-0.34, 0.34))
# first level partition
p1x = [-0.2, 1.0, NaN]; p1y = [1.3, -1.0, NaN]; plot!(p1x, p1y, c = :red, legend = false, width = 3)
# second level partition
p2x = [-1.0, 0.2, NaN, 0.4, 1.2, NaN]; p2y = [-0.8, 0.45, NaN, 0.25, 0.2, NaN];
plot!(p2x, p2y, c = :orange, legend = false, width = 2)
plt = current()
savefig(plt, "../figs/Grid7x3_DAG_VM_NGWP_lvl2_wavelets.png")
