cd(@__DIR__); include("setups/simpletree.jl")
gr(dpi = 200)

## build frames
SGWT = sgwt_frame(Matrix(W); nf = 6)
SGWT = reshape(SGWT, (N, :))
SGWT_dual = (SGWT * SGWT') \ SGWT

distROT = natural_eigdist(𝚽, 𝛌, Q; α = 1.0, input_format = :pmf1, distance = :ROT)
rNGWF, dic_l2x = rngwf_all_vectors(distROT, 𝚽; σ = 0.1 * maximum(distROT), thres = 0.15)
rNGWF_dual = (rNGWF * rNGWF') \ rNGWF
Γ = rngwf_lx(dic_l2x)

## (a)
plt = plot(size = (600, 500))
    gplot!(W, X, width = 1)
    scatter_gplot!(X; marker = f, ms = 4)
    plot!(frame = :none, cbar = true, xlim = [-10, 14])
savefig(plt, "../figs/simpletree_f.png")

## (b)
rel_approx_sgwt, f_approx_sgwt = frame_approx(f, SGWT, SGWT_dual; num_kept = 3 * N)
rel_approx_rngwf, f_approx_rngwf = frame_approx(f, rNGWF, rNGWF_dual; num_kept = 3 * N)
plt = plot(0:length(rel_approx_rngwf)-1, [rel_approx_sgwt rel_approx_rngwf],
    grid = false, lw = 2, c = [:blue :green],
    xlab = "Number of Coefficients Retained",
    ylab = "Relative Approximation Error", yaxis = :log,
    lab = ["SGWT" "rNGWF"])
    plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 11, size = (600, 500))
savefig(plt, "../figs/simpletree_approx_f.png")
