## script for Fig.7, Fig.8(a), Fig.10(a)
using VoronoiDelaunay, VoronoiCells, GeometricalPredicates
using MultiscaleGraphSignalTransforms, Plots, LightGraphs, JLD; gr(dpi = 200)

barbara = JLD.load(joinpath(@__DIR__, "..", "datasets", "barbara_gray_matrix.jld"), "barbara")
G, L, X = SunFlowerGraph(); N = nv(G)
#################### Fig. 7(a) sunflower graph
gplot(1.0*adjacency_matrix(G),X; width=1); scatter_gplot!(X; c = :red, ms = LinRange(1,9,N)); plt = plot!(frame = :none)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/SunFlower.png"))

## Voronoi tessellation
width_x = maximum(abs.(X[:,1])) * 2; width_y = maximum(abs.(X[:,2])) * 2;
width = VoronoiDelaunay.max_coord - VoronoiDelaunay.min_coord
center_coord = (VoronoiDelaunay.min_coord + VoronoiDelaunay.max_coord)/2
X_transform = zeros(N,2)
for i in 1:N
    X_transform[i,:] = X[i,:] ./ [width_x/width, width_y/width] + [center_coord, center_coord]
end

pts = [Point2D(X_transform[i,1], X_transform[i,2]) for i in 1:N]
tess = DelaunayTessellation(N)
push!(tess, pts)
#################### Fig. 7(b) voronoi tessellation
xx, yy = getplotxy(voronoiedges(tess))
plt = plot(xx, yy, xlim=[1,2], ylim=[1,2], linestyle=:auto, linewidth=1, linecolor=:blue, grid=false, label="", aspect_ratio=1, frame=:box)
# savefig(plt, joinpath(@__DIR__, "../paperfigs/Sunflower_Barbara_Voronoi_cells.png"))

#################### Fig. 8(a) barbara sunflower eye
heatmap(barbara, yflip=true, ratio=1, c=:greys); scatter_gplot!(transform2D(X; s = 20, t = [395, 100]); ms = 2, c = :red); sample_location_plt = plot!(cbar = false, frame = :none)
# savefig(sample_location_plt, joinpath(@__DIR__, "../paperfigs/barb_sunflower_eye.png"))

#################### Fig. 10(a) barbara sunflower pants
heatmap(barbara, yflip=true, ratio=1, c=:greys); scatter_gplot!(transform2D(X; s = 20, t = [280, 320]); ms = 2, c = :red); sample_location_plt = plot!(cbar = false, frame = :none)
# savefig(sample_location_plt, joinpath(@__DIR__, "../paperfigs/barb_sunflower_pants.png"))
