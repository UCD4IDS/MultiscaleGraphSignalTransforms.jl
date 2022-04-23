using MultiscaleGraphSignalTransforms, Plots, Graphs, Random
using PrettyTables, StatsPlots

nsim = 500

function generate_ROT_TSD_ratio(nsim, 𝚽, ∇𝚽, 𝛌, Q; edge_length = 1)
    Random.seed!(1234)
    ρ = zeros(nsim)
    for i = 1:nsim
        p = rand(N); p ./= norm(p, 1)
        q = rand(N); q ./= norm(q, 1)
        W1 = ROT_Distance(p, q, Q; edge_length = edge_length)
        K = K_functional(p, q, 𝚽, ∇𝚽, 𝛌; length = edge_length)[1]
        ρ[i] = K / W1
    end
    return ρ
end


function display_basic_stats(ρs)
    header = ["min" "max" "mean" "std"]
    basic_stats = zeros(0, 4)
    for ρ in ρs
        basic_stats = vcat(
            basic_stats,
            round.([minimum(ρ) maximum(ρ) mean(ρ) std(ρ)]; digits = 4)
            )
    end
    pretty_table(basic_stats; header = tuple([header[i] for i = 1:length(header)]))
end



function ROT_TSD_ratio_histogram(ρ)
    plt = histogram(ρ, grid = false, legend = false, c = :teal,
              xlims = [minimum(ρ) - std(ρ), maximum(ρ) + std(ρ)],
              xlab = "ρ")
    return plt
end
