cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

##
frame = sgwt_frame(W; nf = 6)
x = 242
for j = 1:6
    plt = heatmap(reshape(frame[:, x, j], (Nx, Ny))', c = :viridis, ratio = 1,
            frame = :none, xlim = [1, Nx], size = (500, 400))
    savefig(plt, "../figs/Grid$(Nx)x$(Ny)_SGWT_frame_j$(j-1)_x$(x).png")
end
