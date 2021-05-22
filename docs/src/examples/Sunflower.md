# Visualization of the Sunflower Graph

```@example sunflower
using MultiscaleGraphSignalTransforms, LightGraphs, Plots; gr()

# construct the sunflower graph
G, L, X = SunFlowerGraph(); N = nv(G)

# display the sunflower graph (node radii vary for visualization purpose)
gplot(1.0 * adjacency_matrix(G), X; width = 1)
scatter_gplot!(X; c = :red, ms = LinRange(1, 9, N))
plot!(frame = :none, size = (815, 500)) # hide
```
