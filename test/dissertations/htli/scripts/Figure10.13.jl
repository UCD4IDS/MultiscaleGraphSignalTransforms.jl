cd(@__DIR__); include("setups/simpletree.jl")
gr(dpi = 200)

##
ind_eig = [5, 6, 11, 18, 21, 25, 38]
plot(layout = Plots.grid(1, 7), size = (1400, 500))
for i = 1:length(ind_eig)
    l = ind_eig[i]
    w = 𝚽[:, l]
    scatter!(X[:, 1], X[:, 2], marker_z = w, ms = 5, c = :viridis, subplot = i, mswidth = 0)
    plot!(frame = :none, cbar = false, grid = false, legend = false, subplot = i)
    plot!(title = latexstring("\\mathbf{\\phi}_{", l - 1, "}"), titlefont = 18, subplot = i)
end
plot!(subplot = 10, grid = false, frame = :none)
plt = current()
savefig(plt, "../figs/simpletree_f_rngwf_top21_ref_eigenvectors.png")
