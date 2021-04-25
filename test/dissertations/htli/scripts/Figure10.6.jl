cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

##
x = 242

## (a)
l = 11 # vertical
plt = grid_eigenvector_plot(l)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_eigenvector_l$(l-1).png")

## (d)
plt = grid_NGWFvector_plot(l, x, D, 𝚽)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_NGWF_DAG_sig02_l$(l-1)_x$(x).png")

## (b)
l = 16  # horizontal
plt = grid_eigenvector_plot(l)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_eigenvector_l$(l-1).png")

## (e)
plt = grid_NGWFvector_plot(l, x, D, 𝚽)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_NGWF_DAG_sig02_l$(l-1)_x$(x).png")

## (c)
l = 30  # mix
plt = grid_eigenvector_plot(l)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_eigenvector_l$(l-1).png")

## (f)
plt = grid_NGWFvector_plot(l, x, D, 𝚽)
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_NGWF_DAG_sig02_l$(l-1)_x$(x).png")
