cd(@__DIR__); include("setups/simpletree.jl")
gr(dpi = 200)

## build rNGWF
distROT = natural_eigdist(𝚽, 𝛌, Q; α = 1.0, input_format = :pmf1, distance = :ROT)
rNGWF, dic_l2x = rngwf_all_vectors(distROT, 𝚽; σ = 0.1 * maximum(distROT), thres = 0.15)
Γ = rngwf_lx(dic_l2x)

just_eigenvec = [2, 16]
important_idx = sortperm((rNGWF' * f).^2; rev = true)
plot(layout = Plots.grid(3, 7), size = (1400, 1500))
for i = 1:21
    l, x = Γ[important_idx[i]]
    w = rNGWF[:, important_idx[i]]
    w ./= norm(w, 2)
    scatter!(X[:, 1], X[:, 2], marker_z = w, ms = 5, c = :viridis, subplot = i, mswidth = 0)
    plot!(frame = :none, cbar = false, grid = false, legend = false, subplot = i)
    if l in just_eigenvec
        plot!(title = latexstring("\\bar{\\psi}_{", l, ",", x, "} \\equiv \\mathbf{\\phi}_{$(l)}"), titlefont = 18, subplot = i)
    else
        plot!(title = latexstring("\\bar{\\psi}_{", l, ",", x, "}"), titlefont = 18, subplot = i)
    end
end
plt = current()
savefig(plt, "../figs/simpletree_f_rngwf_top21.png")
