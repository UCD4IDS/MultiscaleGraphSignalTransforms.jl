cd(@__DIR__); include("setups/rgc100.jl")
gr(dpi = 200)

plt = plot(0:(N - 1), 𝛌, c = :black, lw = 1, ylims = [0, 5], legend = false,
        shape = :circle, frame = :box, xlab = latexstring("l"),
        ylab = latexstring("\\tilde{\\lambda}_l"), size = (500, 400))
savefig(plt, "../figs/RGC100_eigenvalues.png")
