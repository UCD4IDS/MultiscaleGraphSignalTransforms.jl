cd(@__DIR__); include("setups/grid7x3.jl")
pyplot(dpi = 200)

## (a) eigenvectors by nondecreasing eigenvalue ordering
plot(layout = Plots.grid(3, 7))
for i in 1:N
    heatmap!(reshape(𝚽[:, i], (Nx, Ny))', c = :viridis, cbar = false,
                clims = (-0.4,0.4), frame = :none, ratio = 1, ylim = [0, Ny + 1],
                title = latexstring("\\phi_{", i-1, "}"), titlefont = 12,
                subplot = i)
end
plt = current()
savefig(plt, "../figs/grid7x3_evsp_title.png")

## (b) eigenvectors by natural frequency ordering
plot(layout = Plots.grid(3, 7))
for i in 1:N
    k = grid2eig_ind[i]
    heatmap!(reshape(𝚽[:,k], (Nx, Ny))', c = :viridis, cbar = false,
                clims = (-0.4,0.4), frame = :none, ratio = 1, ylim = [0, Ny + 1],
                title = latexstring("\\varphi_{", string(eig2dct[k,1]),
                ",", string(eig2dct[k,2]), "}"), titlefont = 12, subplot = i)
end
plt = current()
savefig(plt, "../figs/grid7x3_dct_title2.png")
