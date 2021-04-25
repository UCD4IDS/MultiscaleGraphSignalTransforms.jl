cd(@__DIR__); include("setups/path512.jl")
pyplot(dpi = 200)

##
@time VM_NGWP = vm_ngwp(𝚽, GP_dual)
@time LP_NGWP = lp_ngwp(𝚽, Gstar_Sig.W, GP_dual; ϵ = 0.3)

dict_VM = Dict()
dict_LP = Dict()
for j in [3, 4, 5]
   for k in [1, 2, 3]
      if GP_dual.rs[k + 1, j] - GP_dual.rs[k, j] < 1
         continue
      end
      WW_VM = sort_wavelets(VM_NGWP[GP_dual.rs[k, j]:(GP_dual.rs[k + 1, j] - 1), j, :]')
      WW_LP = sort_wavelets(LP_NGWP[GP_dual.rs[k, j]:(GP_dual.rs[k + 1, j] - 1), j, :]')
      Y = [-0.5, 0.6]
      l = Int(floor(size(WW_LP, 2) / 2))
      dict_VM[(j-1, k-1, l-1)] = WW_VM[:, l]
      dict_LP[(j-1, k-1, l-1)] = WW_LP[:, l]
      plt = plot(size = (400, 350), layout = Plots.grid(2, 1), ylim = Y)
         plot!(WW_VM[:, l], c = :black, grid = false, frame = :box, lw = 0.8, xtick = false,
         lab = latexstring("\\psi^{($(j-1))}_{$(k-1), $(l-1)} (\\mathrm{VM})"), legendfontsize = 11)
         plot!(WW_LP[:, l], c = :black, grid = false, frame = :box, lw = 0.8,
         lab = latexstring("\\psi^{($(j-1))}_{$(k-1), $(l-1)} (\\mathrm{LP})"), legendfontsize = 11,
         subplot = 2)
         xticks!([1; 64:64:N], vcat(string(1), [string(k) for k in 64:64:N]),
         subplot = 2)
      savefig(plt, "../figs/Path512_VM_LP_j$(j-1)k$(k-1)l$(l-1).png")
   end
end

## Table 8.2
using PrettyTables
header = ["basis vector" "main support width" "sidelobe attenuation"]
res = Matrix{String}(undef, 0, 3)
iter = sort(collect(keys(dict_VM)))
for γ in iter
   ms = find_mainsupport(dict_VM[γ]; ϵ = 0.01)
   sa = sidelobe_attenuation(dict_VM[γ])
   j, k, l = γ
   res = vcat(res, ["ψ^{$(j)}_{$(k), $(l)} (VM)" ms[2] - ms[1] + 1 round(sa; digits = 4)])
   ms = find_mainsupport(dict_LP[γ]; ϵ = 0.01)
   sa = sidelobe_attenuation(dict_LP[γ])
   res = vcat(res, ["ψ^{$(j)}_{$(k), $(l)} (LP)" ms[2] - ms[1] + 1 round(sa; digits = 4)])
end
pretty_table(res, header; alignment = :c)
