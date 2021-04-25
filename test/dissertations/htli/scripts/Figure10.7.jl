cd(@__DIR__); include("setups/grid23x22.jl")
gr(dpi = 200)

## assemble frames
frame = sgwt_frame(W; nf = 6)
SGWT = reshape(frame, (N, :))
SGWT_dual = (SGWT * SGWT') \ SGWT

rNGWF, dic_l2x = rngwf_all_vectors(D, 𝚽; σ = 0.2 * maximum(D), thres = 0.15)
rNGWF_dual = (rNGWF * rNGWF') \ rNGWF
Γ = rngwf_lx(dic_l2x)

## (a)
f = digit_img[:]
plt = heatmap(reshape(f, (Nx, Ny))', c = :viridis, ratio = 1,
        frame = :none, xlim = [1, Nx], size = (600, 500))
savefig(plt, "../figs/grid23x22_fdigit3.png")

## (b)
rel_approx_sgwt, f_approx_sgwt = frame_approx(f, SGWT, SGWT_dual; num_kept = 3 * N)
rel_approx_rngwf, f_approx_rngwf = frame_approx(f, rNGWF, rNGWF_dual; num_kept = 3 * N)
plot(0:length(rel_approx_rngwf)-1, [rel_approx_sgwt rel_approx_rngwf],
    grid = false, lw = 2, c = [:blue :green],
    xlab = "Number of Coefficients Retained",
    ylab = "Relative Approximation Error", yaxis = :log,
    lab = ["SGWT" "rNGWF"])
plt = plot!(xguidefontsize = 14, yguidefontsize = 14, legendfontsize = 11, size = (600, 500))
savefig(plt, "../figs/grid23x22_approx_fdigit3.png")
