# Visualization of the Sunflower Graph
```@docs
SunFlowerGraph
```
Let us see how to visualize the sunflower graph by `gplot()` and `scatter_gplot()`.
```@example sunflower
using MultiscaleGraphSignalTransforms, LightGraphs, Plots

# construct the sunflower graph
G, L, X = SunFlowerGraph(); N = nv(G)

# display the sunflower graph (node radii vary for visualization purpose)
gplot(1.0 * adjacency_matrix(G), X; width = 1)
scatter_gplot!(X; c = :red, ms = LinRange(1, 9, N))
plot!(frame = :none, size = (815, 500)) # hide
```
One can also represent a signal on the graph by colors. For example,
```@example sunflower
f = zeros(N)
f[1:200] .= 1
f[301:N] .= -1

# display the graph signal
gplot(1.0 * adjacency_matrix(G), X; width = 1)
scatter_gplot!(X; marker = f, ms = LinRange(1, 9, N))
plot!(frame = :none, size = (815, 500)) # hide
```
