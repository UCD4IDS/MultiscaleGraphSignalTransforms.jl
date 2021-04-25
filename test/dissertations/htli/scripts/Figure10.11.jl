cd(@__DIR__); include("setups/simpletree.jl")
gr(dpi = 200)

## build SGWT frame
SGWT = sgwt_frame(Matrix(W); nf = 6)
SGWT = reshape(SGWT, (N, :))

important_idx = sortperm((SGWT' * f).^2; rev = true)
plot(layout = Plots.grid(3, 7), size = (1400, 1500))
for i = 1:21
    ind = important_idx[i]
    j = Int(floor(ind / N))
    x = ind - j * N
    w = SGWT[:, ind]
    w ./= norm(w, 2)
    scatter!(X[:, 1], X[:, 2], marker_z = w, ms = 5, c = :viridis, subplot = i, mswidth = 0)
    plot!(frame = :none, cbar = false, grid = false, legend = false, subplot = i)
    plot!(title = latexstring("\\psi^{\\mathrm{SGWT}}_{", j, ",", x, "}"), titlefont = 18, subplot = i)
end
plt = current()
savefig(plt, "../figs/simpletree_f_sgwt_top21.png")
