cd(@__DIR__); include("setups/path128.jl")
gr(dpi = 200)

l = 64
𝛍 = nat_spec_filter(l, distDCT; σ = 0.05 * maximum(distDCT))
plot(0:N-1, 𝛍, legend = true, c = :blue,
     grid = false, frame = :box, marker = :circle, ms = 3, lw = 1, mswidth = 0,
     lab = latexstring("\\sigma=0.05 \\cdot d_{\\max}"))

𝛍 = nat_spec_filter(l, distDCT; σ = 0.1 * maximum(distDCT))
plot!(0:N-1, 𝛍, legend = true, c = :red,
     grid = false, frame = :box, marker = :circle, ms = 3, lw = 1, mswidth = 0,
     lab = latexstring("\\sigma=0.1 \\cdot d_{\\max}"))

𝛍 = nat_spec_filter(l, distDCT; σ = 0.2 * maximum(distDCT))
plot!(0:N-1, 𝛍, ylim = [0, 0.1], legend = true, c = :green,
     grid = false, frame = :box, marker = :circle, ms = 3, lw = 1, mswidth = 0,
     lab = latexstring("\\sigma=0.2 \\cdot d_{\\max}"), tickfontsize=10, legendfontsize = 12)
xticks!([0; 15:16:127], vcat(string("DC"), [string(k) for k in 15:16:127]))
plot!(xlab = latexstring("i"), ylab = latexstring("\\mu_{63}"), guidefontsize = 16)
plt = current()
savefig(plt, "../figs/Path128_NSGfilters_mu$(l-1).png")
display(plt)
