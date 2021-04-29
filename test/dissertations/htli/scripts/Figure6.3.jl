cd(@__DIR__); include("setups/path512.jl");
gr(dpi = 200)

## Build VM-NGWP
@time VM_NGWP = vm_ngwp(𝚽, GP_dual)

#################### Fig.5
j = 5
for k in [1, 2, 5]
   WW = sort_wavelets(VM_NGWP[GP_dual.rs[k, j]:(GP_dual.rs[k + 1, j] - 1), j, :]')
   sc = 0.75
   if k == 1
      sc = 0.5
   end
   if k == 2
      WW[:, end] *= -1
   end
   plt = wiggle(WW; sc = sc)
   savefig(plt, "../figs/Path512_VM_NGWP_j$(j-1)k$(k-1).png")
end
current()
