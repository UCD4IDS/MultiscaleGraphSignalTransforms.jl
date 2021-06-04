using Plots, SparseArrays, JLD2, LinearAlgebra, MultiscaleGraphSignalTransforms

#@load "MN_MutGauss.jld2"
#tmp1 = G["G"]
#G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])
#color_limit = (-0.2,2.2)

JLD2.@load "Toronto.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])

GraphSig_Plot(G, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0.)



G = Adj2InvEuc(G)
GP = partition_tree_fiedler(G,:Lrw)
dmatrix = ghwt_analysis!(G, GP=GP)


############# Haar
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

############# Walsh
BS_walsh = bs_walsh(GP)
dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)

############# GHWT_c2f
dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)

############# GHWT_f2c
dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)

############# eGHWT
dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)



################################################################################
####################### Approximation error plot################################
################################################################################
function approx_error(DVEC::Array{Array{Float64,1},1})
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    frac = 0:0.01:0.3
    T = ["Haar","Walsh","GHWT_c2f", "GHWT_f2c", "eGHWT"]
    L = [(:dashdot,:orange),(:dashdot,:blue),(:solid, :red),(:solid, :green),(:solid, :black)]
    for i = 1:5
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2, rev = true)

        er = fill(0., length(frac))

        for j = 1:length(frac)
            p = Int64(floor(frac[j]*N))
            er[j] = sqrt(dvec_norm^2 - sum(dvec_sort[1:p]))/dvec_norm
        end

        plot!(frac, er, yaxis=:log, xlims = (0.,0.3), label = T[i], line = L[i])
        print(er[26])
    end
end




################################################################################
######################### Synthesized by top  vectors###########################
################################################################################
color_limit_residual = (0., 0.15)

function top_vectors_residual(p::Int64, dvec::Array{Float64,2}, BS::BasisSpec, GP::GraphPart, G::GraphSig)
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0
    (f, GS) = ghwt_synthesis(dvecT, GP, BS, G)
    GS.name = "Square of the residual"
    #GS.f =  (G.f - GS.f).^2
    #print((maximum(GS.f), minimum(GS.f)))
    #GraphSig_Plot(GS)
    GS.f = (G.f - GS.f).^2 ./(G.f.^2)
    GraphSig_Plot(GS, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., clim = color_limit_residual)
    current()
end





################################
### Generate results
################################
# Figure 3(b)
approx_error([dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:]])
current()

###
frac = 1/4
p = Int64(ceil(frac*G.length))

### original (Figure )
# Figure 3(a)
GraphSig_Plot(G, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., clim = (2000., 110000.))

### haar
# Figure 4(a)
top_vectors_residual(p, dvec_haar, BS_haar, GP, G)

### walsh
# Figure 4(b)
top_vectors_residual(p, dvec_walsh, BS_walsh, GP, G)

### c2f
# Figure 4(b)
top_vectors_residual(p, dvec_c2f, BS_c2f, GP, G)

### f2c
# Figure 4(c)
top_vectors_residual(p, dvec_f2c, BS_f2c, GP, G)

### eghwt
# Figure 4(d)
top_vectors_residual(p, dvec_eghwt, BS_eghwt, GP, G)
