cd(@__DIR__); allNGWPs = false; include("setups/toronto.jl")
gr(dpi = 200)

## fdensity 16 most important VM-NGWP vectors (ignore the DC vector)
VM_NGWP = vm_ngwp(𝚽, GP_dual)
f = zeros(N); for i in 1:N; f[i] = length(findall(dist_X[:,i] .< 1/minimum(edge_weight))); end
G_Sig.f = reshape(f, (N, 1))
dmatrix_VM = ngwp_analysis(G_Sig, VM_NGWP)
dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
important_idx = sortperm(dvec_vm_ngwp[:].^2; rev = true)
println("================fdensity-VM-NGWP-top-basis-vectors=================")
for i in 2:17
    dr, dc = BS_vm_ngwp.levlist[important_idx[i]]
    w = VM_NGWP[dr, dc, :]
    j, k, l = NGWP_jkl(GP_dual, dr, dc)
    print("(j, k, l) = ($(j), $(k), $(l))  ")
    if j == jmax
        print("φ_{$(GP_dual.ind[dr]-1)}")
    end
    println()
    gplot(A, X; width=1)
    scatter_gplot!(X; marker = w, plotOrder = :s2l, ms = 3)
    plt = plot!(cbar = false, clims = (-0.075,0.075))
    savefig(plt,
        "../figs/Toronto_fdensity_DAG_VM_NGWP_ibv$(lpad(i,2,"0"))_j$(j)_k$(k)_l$(l).png")
end
