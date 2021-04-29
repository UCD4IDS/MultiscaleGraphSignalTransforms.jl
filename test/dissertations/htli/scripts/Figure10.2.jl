cd(@__DIR__); include("setups/path128.jl")
gr(dpi = 200)

##
f = 1.5 * 𝚽[:, 11] + 𝚽[:, 31] .* characteristic(1:32, N) +
    𝚽[:, 61] .* characteristic(1:64, N) + 0.5 * 𝚽[:, 111]

## (a)
plot(f, c = :black, grid = false, legend = false, frame = :box, lw = 1.5)
xticks!([1; 16:16:128], vcat(string(1), [string(k) for k in 16:16:128]))
plt = plot!(xlab = latexstring("x"), ylab = latexstring("f"),
    size = (600, 400), guidefontsize = 18, xlim = [1, N])
savefig(plt, "../figs/Path128_Spectrogram_f.png")

## (b)
plot(0:N-1, 𝚽' * f, c = :black, grid = false, legend = false, frame = :box, lw = 1.5)
xticks!([0; 15:16:127], vcat(string("DC"), [string(k) for k in 15:16:127]))
plt = plot!(xlab = latexstring("l"), ylab = latexstring("g"),
    size = (600, 400), guidefontsize = 18, xlim = [1, N])
savefig(plt, "../figs/Path128_Spectrogram_g.png")
