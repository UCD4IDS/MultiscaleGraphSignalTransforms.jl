cd(@__DIR__); include("setups/path128.jl")
gr(dpi = 200)

##
f = 1.5 * 𝚽[:, 11] + 𝚽[:, 31] .* characteristic(1:32, N) +
    𝚽[:, 61] .* characteristic(1:64, N) + 0.5 * 𝚽[:, 111]

## (a)
plt = path_spectrogram(f, distDCT, 𝚽; c = 0.01)
savefig(plt, "../figs/Path128_Spectrogram_CoeffEnergy_sig001dmax.png")

## (b)
plt = path_spectrogram(f, distDCT, 𝚽; c = 0.02)
savefig(plt, "../figs/Path128_Spectrogram_CoeffEnergy_sig002dmax.png")

## (c)
plt = path_spectrogram(f, distDCT, 𝚽; c = 0.05)
savefig(plt, "../figs/Path128_Spectrogram_CoeffEnergy_sig005dmax.png")

## (d)
plt = path_spectrogram(f, distDCT, 𝚽; c = 0.1)
savefig(plt, "../figs/Path128_Spectrogram_CoeffEnergy_sig01dmax.png")
