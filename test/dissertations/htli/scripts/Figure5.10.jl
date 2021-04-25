cd(@__DIR__); include("setups/simpletree.jl");
gr(dpi=200)

plt = plot(0:(N - 1), 𝛌, c = :black, lw = 1, ylims = [0, 5], legend = false,
        shape = :circle, frame = :box, xlab = latexstring("l"), ms = 4,
        ylab = latexstring("\\lambda_l"), size = (500, 400))
savefig(plt, "../figs/SimpleTree_eigenvalues.png")
