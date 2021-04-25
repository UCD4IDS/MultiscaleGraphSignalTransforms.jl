cd(@__DIR__); runapprox = true; allNGWPs = false; include("setups/sunflower.jl")
gr(dpi = 200)

##
VM_NGWP_Lsym = vm_ngwp(𝚽sym, GP_dual_Lsym)
f = matread("../datasets/sunflower_barbara_voronoi.mat")["f_eye_voronoi"]
G_Sig.f = reshape(f, (N, 1))
dmatrix_VM_Lsym = ngwp_analysis(G_Sig, VM_NGWP_Lsym)
dvec_vm_ngwp_Lsym, BS_vm_ngwp_Lsym = ngwp_bestbasis(dmatrix_VM_Lsym, GP_dual_Lsym)
important_idx = sortperm(dvec_vm_ngwp_Lsym[:].^2; rev = true)
println("================feye-VM-NGWP-Lsym-top-basis-vectors=================")
for i in 2:17
    dr, dc = BS_vm_ngwp_Lsym.levlist[important_idx[i]]
    w = VM_NGWP_Lsym[dr, dc, :]
    j, k, l = NGWP_jkl(GP_dual_Lsym, dr, dc)
    print("(j, k, l) = ($(j), $(k), $(l))  ")
    if j == jmax_Lsym
        print("φ_{$(GP_dual_Lsym.ind[dr]-1)}")
    end
    println()
    scatter_gplot(X; marker = w, ms = LinRange(4.0, 14.0, N), c = :greys)
    plt = plot!(xlim = [-1.2, 1.2], ylim = [-1.2, 1.2], frame = :none,
                cbar = false, clims = (-0.15, 0.15))
    savefig(plt,
        "../figs/SunFlower_feye_DAG_VM_NGWP_Lsym_ibv$(lpad(i,2,"0"))_j$(j)_k$(k)_l$(l).png")
end
