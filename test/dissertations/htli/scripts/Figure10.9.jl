cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

## build rNGWF
rNGWF, dic_l2x = rngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D), thres = 0.15)
rNGWF_dual = (rNGWF * rNGWF') \ rNGWF
Γ = rngwf_lx(dic_l2x)
f = digit_img[:]

important_idx = sortperm((rNGWF' * f).^2; rev = true)
red_box_inds = [41, 43, 66, 69, 71, 73, 75, 84, 85, 86, 89, 99]
plot(layout = Plots.grid(10,10), size = (2300, 2200))
for i = 1:100
    plot!(clims = (-0.002, 0.02))
    grid_vector_plot!(important_idx[i], i, rNGWF)
    if i in red_box_inds
        plot_square!(Nx, Ny; subplot = i)
    end
end
plt = current()
savefig(plt, "../figs/Grid23x22_fdigit3_rngwf_top100.png")
