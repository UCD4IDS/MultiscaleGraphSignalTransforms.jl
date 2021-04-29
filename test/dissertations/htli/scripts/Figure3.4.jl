cd(@__DIR__); include("setups/rgc100.jl")
gr(dpi = 200)

plt = gplot(W, X; width = 1); scatter_gplot!(X; marker = 𝚽[:, 1142], ms = 4)
        plot!(xlims = [-180, 220], xlabel = "x(μm)", ylabel = "y(μm)",
        clims = (-0.12, 0.12), frame = :box, size = (550, 500), right_margin = 5mm)
savefig(plt, "../figs/RGC100_twotypes_semi_oscillation.png")

plt = gplot(W, X; width = 1); scatter_gplot!(X; marker = 𝚽[:, 1143], ms = 4,
        plotOrder = :s2l); plot!(xlims = [-180, 220], xlabel = "x(μm)",
        ylabel = "y(μm)", clims = (-0.3, 0.3), frame = :box, size = (550, 500),
        right_margin = 5mm)
savefig(plt, "../figs/RGC100_twotypes_localized.png")
