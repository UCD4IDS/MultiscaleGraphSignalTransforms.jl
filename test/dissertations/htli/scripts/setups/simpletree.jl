using MultiscaleGraphSignalTransforms, JLD, Plots, LightGraphs, MultivariateStats
using LaTeXStrings
push!(LOAD_PATH, @__DIR__)
using pSGWT

G = loadgraph("../datasets/simple_tree_graph.lgz")
X = load("../datasets/simple_tree_xy.jld", "xy")
N = nv(G)
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L)
standardize_eigenvectors!(𝚽)
W = 1.0 * adjacency_matrix(G)
Q = incidence_matrix(G; oriented = true)

ib1 = 36:56
ib2 = 21:35
ib3 = 71:100
ib4 = 57:70
ijc = [3,5,12,16]
ir = setdiff(1:N, ib1, ib2, ib3, ib4, ijc)


##
P = 𝚽.^2
# P = exp.(𝚽) ./ sum(exp.(𝚽), dims = 1)

dist_sROT, Ws, Xs, 𝚯 = eigsROT_Distance(P, W, X; α = 1.0)

##
function find_active_eigenvectors(P, interest_locs; threshold = 0.5)
    N = size(P, 1)
    energy = zeros(N)
    for k=1:N
        energy[k] = norm(P[interest_locs, k], 1) / norm(P[:, k], 1)
    end
    ind = findall(energy .> threshold)
    return ind
end

# index of eigenvectors active at branch k (k = 1,2,3,4)
ieb1 = find_active_eigenvectors(𝚽.^2, ib1)
ieb2 = find_active_eigenvectors(𝚽.^2, ib2)
ieb3 = find_active_eigenvectors(𝚽.^2, ib3)
ieb4 = find_active_eigenvectors(𝚽.^2, ib4)
iejc = find_active_eigenvectors(𝚽.^2, ijc; threshold = 0.1)


function simpletree_mds_plot(E, ieb1, ieb2, ieb3, ieb4, iejc)
    scatter_gplot(E'; c = :grey, ms = 2)
    scatter_gplot!(E[:, ieb1]'; c = :pink, ms = 2)
    scatter_gplot!(E[:, ieb2]'; c = :orange, ms = 2)
    scatter_gplot!(E[:, ieb3]'; c = :green, ms = 2)
    scatter_gplot!(E[:, ieb4]'; c = :yellow, ms = 2)
    scatter_gplot!(E[:, iejc]'; c = :red, ms = 2)
    scatter_gplot!(E[:, 1:1]'; c = :magenta, ms = 4)
    plt = plot!(xaxis = "X₁", yaxis = "X₂", zaxis = "X₃", legend = false,
                cbar = false, grid = true)
    return plt
end

##
f = zeros(N)
ib1 = 36:56; ib2 = 21:35; ib3 = 71:100; ib4 = 57:70; ijc = [3,5,12,16]
f[ib1] = sin.(0.3 * (0:(length(ib1) - 1)))
f[ib2] = cos.(0.4 * (0:(length(ib2) - 1)))
f[ib3] = sin.(0.5 * (0:(length(ib3) - 1)))
f[ib4] = cos.(0.6 * (0:(length(ib4) - 1)))
f[ijc] .= 1;
