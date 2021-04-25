cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

## build SGWT frame
frame = sgwt_frame(W; nf = 6)
SGWT = reshape(frame, (N, :))
SGWT_dual = (SGWT * SGWT') \ SGWT
f = digit_img[:]

important_idx = sortperm((SGWT' * f).^2; rev = true)
plot(layout = Plots.grid(10,10), size = (2300, 2200))
for i = 1:100
    grid_vector_plot!(important_idx[i], i, SGWT)
    if i == 90
        plot_square!(Nx, Ny; subplot = i)
    end
end
plt = current()
savefig(plt, "../figs/Grid23x22_fdigit3_sgwt_top100.png")
