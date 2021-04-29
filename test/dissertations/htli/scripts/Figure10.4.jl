cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

##
frame = sgwt_frame(W; nf = 6)
filters = sgwt_filter_banks(W, 𝛌; nf = 6)
plot(𝛌, filters', lw = 3, size = (800, 450), grid = false, frame = :box,
    legendfontsize = 14, xlab = latexstring("\\lambda"), xguidefontsize = 14,
    lab = [latexstring("h(\\lambda)") latexstring("g(s_{1}\\lambda)") latexstring("g(s_{2}\\lambda)") latexstring("g(s_{3}\\lambda)") latexstring("g(s_{4}\\lambda)") latexstring("g(s_{5}\\lambda)")])
# plot!(𝛌, sum(filters, dims = 1)[:]/6, c = :black, lw = 3, lab = "partition of unity")
plt = current()
savefig(plt, "../figs/Grid$(Nx)x$(Ny)_SGWT_filter_banks.png")
