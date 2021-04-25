cd(@__DIR__); runapprox = true; allNGWPs = false; include("setups/sunflower.jl")
gr(dpi = 200)

##
VM_NGWP = vm_ngwp(𝚽, GP_dual)
f = matread("../datasets/sunflower_barbara_voronoi.mat")["f_trouser_voronoi"]
G_Sig.f = reshape(f, (N, 1))
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
println("================ftrouser-VM-NGWP-top-basis-vectors=================")
for i in 2:17
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    j, k, l = NGWP_jkl(GP_dual, dr, dc)
    print("(j, k, l) = ($(j), $(k), $(l))  ")
    if j == jmax
        print("φ_{$(GP_dual.ind[dr]-1)}")
    end
    println()
    scatter_gplot(X; marker = w, ms = LinRange(4.0, 14.0, N), c = :greys)
    plt = plot!(xlim = [-1.2, 1.2], ylim = [-1.2, 1.2], frame = :none,
                cbar = false, clims = (-0.15, 0.15))
    savefig(plt,
        "../figs/SunFlower_ftrouser_DAG_VM_NGWP_ibv$(lpad(i,2,"0"))_j$(j)_k$(k)_l$(l).png")
end
