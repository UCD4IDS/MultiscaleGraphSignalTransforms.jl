using NGWP, Plots, LightGraphs, JLD, MTSG, MAT

barbara = JLD.load("../datasets/barbara_gray_matrix.jld", "barbara")

## Build weighted graph
G, L, X = SunFlowerGraph(N = 400); N = nv(G)
W = 1.0 * adjacency_matrix(G)
if runapprox
    Q = incidence_matrix(G; oriented = true)
    edge_weight = [e.weight for e in edges(G)]

    ## eigenvectors of L(G)
    𝛌, 𝚽 = eigen(Matrix(L))
    standardize_eigenvectors!(𝚽)

    ## eigenvectors of Lsym(G)
    deg = sum(W, dims = 1)[:]  # weighted degree vector
    Lsym = diagm(deg.^(-1/2)) * (diagm(deg) - W) * diagm(deg.^(-1/2))
    𝛌sym, 𝚽sym = eigen(Lsym)
    standardize_eigenvectors!(𝚽sym)

    ## Build Dual Graph by DAG metric
    distDAG = eigDAG_Distance(𝚽, Q, N; edge_weight = edge_weight)
    Gstar_Sig = dualgraph(distDAG)
    G_Sig = GraphSig(W, xy = X)
    GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
    GP_primal = pairclustering(𝚽, GP_dual)
    jmax = size(GP_dual.rs, 2) - 1  # zero-indexed

    if allNGWPs
        # 54.986524 seconds (1.19 M allocations: 25.127 GiB, 2.36% gc time)
        @time VM_NGWP = vm_ngwp(𝚽, GP_dual)
        # 0.611176 seconds (225.41 k allocations: 844.488 MiB, 11.71% gc time)
        @time PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal)
        # 50.939912 seconds (7.67 M allocations: 28.051 GiB, 3.11% gc time)
        @time LP_NGWP = lp_ngwp(𝚽, Gstar_Sig.W, GP_dual; ϵ = 0.3)
    end

    ## Build Dual Graph by DAG metric (Lsym)
    distDAG_Lsym = eigDAG_Distance(𝚽sym, Q, N; edge_weight = edge_weight)
    Gstar_Sig_Lsym = dualgraph(distDAG_Lsym)
    GP_dual_Lsym = partition_tree_fiedler(Gstar_Sig_Lsym; swapRegion = false)
    GP_primal_Lsym = pairclustering(𝚽sym, GP_dual_Lsym)
    jmax_Lsym = size(GP_dual_Lsym.rs, 2) - 1

    if allNGWPs
        # 61.328441 seconds (1.34 M allocations: 29.058 GiB, 2.30% gc time)
        @time VM_NGWP_Lsym = vm_ngwp(𝚽sym, GP_dual_Lsym)
        # 0.605507 seconds (237.87 k allocations: 876.689 MiB, 10.83% gc time)
        @time PC_NGWP_Lsym = pc_ngwp(𝚽sym, GP_dual_Lsym, GP_primal_Lsym)
        # 77.918189 seconds (2.12 M allocations: 42.872 GiB, 2.91% gc time)
        @time LP_NGWP_Lsym = lp_ngwp(𝚽sym, Gstar_Sig_Lsym.W, GP_dual_Lsym; ϵ = 0.3)
    end
else
    using VoronoiDelaunay, VoronoiCells, GeometricalPredicates
    ## Voronoi tessellation
    width_x = maximum(abs.(X[:, 1])) * 2; width_y = maximum(abs.(X[:, 2])) * 2;
    width = VoronoiDelaunay.max_coord - VoronoiDelaunay.min_coord
    center_coord = (VoronoiDelaunay.min_coord + VoronoiDelaunay.max_coord)/2
    X_transform = zeros(N,2)
    for i in 1:N
        X_transform[i,:] = X[i,:] ./ [width_x/width, width_y/width] + [center_coord, center_coord]
    end

    pts = [Point2D(X_transform[i,1], X_transform[i,2]) for i in 1:N]
    tess = DelaunayTessellation(N)
    push!(tess, pts)
end
