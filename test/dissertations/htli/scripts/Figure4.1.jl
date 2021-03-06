cd(@__DIR__); include("setups/distROT_distTSD_ratio.jl")
gr(dpi = 200)

## (a) P64
include("setups/path64.jl")
Ï1 = generate_ROT_TSD_ratio(nsim, ð½, âð½, ð, Q)
# plt = ROT_TSD_ratio_histogram(Ï1)
# savefig(plt, "../figs/Path64_ROT_TSD.png")

## (b) Pâ x Pâ
include("setups/grid7x3.jl")
Ï2 = generate_ROT_TSD_ratio(nsim, ð½, âð½, ð, Q)
# plt = ROT_TSD_ratio_histogram(Ï2)
# savefig(plt, "../figs/Grid7x3_ROT_TSD.png")

## (c) Erdos RÃ©nyi
include("setups/er.jl")
Ï3 = generate_ROT_TSD_ratio(nsim, ð½, âð½, ð, Q)
# plt = ROT_TSD_ratio_histogram(Ï3)
# savefig(plt, "../figs/Erdos_Renyi_ROT_TSD.png")

## (d) weighted RGC100
include("setups/rgc100.jl")
Ï4 = generate_ROT_TSD_ratio(nsim, ð½, âð½, ð, Q; edge_length = edge_length)
# plt = ROT_TSD_ratio_histogram(Ï4)
# savefig(plt, "../figs/wRGC100_ROT_TSD.png")

## boxplot
plt = boxplot(["(a) Path" "(b) Grid" "(c) ER" "(d) RGC100"], [Ï1, Ï2, Ï3, Ï4];
    legend = false, frame = :box, ylim = [0.9, 1.9], tickfontsize = 11,
    outliers = true, grid = false, range = 3, lw = 1, size = (800, 600),
    ylab = "Ï", yguidefontsize = 14, xtickfontsize = 12, color = :white)
savefig(plt, "../figs/ROT_TSD_boxplots.png")

## Table 4.1
display_basic_stats([Ï1, Ï2, Ï3, Ï4])
