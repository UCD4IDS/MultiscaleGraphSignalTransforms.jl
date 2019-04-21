using FileIO, Plots, SparseArrays, JLD2, LinearAlgebra, MTSG

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
    end
end

approx_error([dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:]])
current()

#Plots.savefig("approx_error.pdf")







################################################################################
######################### Synthesized by top  vectors###########################
################################################################################
### Plot the graph synthesized by top frac*total_number vectors
color_limit = (2000., 110000.)
function top_vectors_synthesis(p::Int64, dvec::Array{Float64,2}, BS::BasisSpec, GP::GraphPart, G::GraphSig)
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0
    (f, GS) = ghwt_synthesis(dvecT, GP, BS, G)
    #print((maximum(f), minimum(f)))
    GS.name = "Synthesized by top "*string(p)*" vectors"
    #GraphSig_Plot(GS)
    GraphSig_Plot(GS, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., clim = color_limit)
    current()
end

### Plot the residual squared
color_limit_residual = (0., 1e8)
function top_vectors_residual(p::Int64, dvec::Array{Float64,2}, BS::BasisSpec, GP::GraphPart, G::GraphSig)
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0
    (f, GS) = ghwt_synthesis(dvecT, GP, BS, G)
    GS.name = "Square of the residual"
    GS.f =  (G.f - GS.f).^2
    #print((maximum(GS.f), minimum(GS.f)))
    #GraphSig_Plot(GS)
    GraphSig_Plot(GS, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., clim = color_limit_residual)
    current()
end

### Plot the histogram of the residual squared
function top_vectors_residual_hist(p::Int64, dvec::Array{Float64,2}, BS::BasisSpec, GP::GraphPart, G::GraphSig)
    sorted_dvec = sort(abs.(dvec[:]), rev = true)
    dvecT = copy(dvec)
    dvecT[abs.(dvec) .< sorted_dvec[p]].= 0
    (f, GS) = ghwt_synthesis(dvecT, GP, BS, G)
    GS.name = "Residual of the approximation signal"
    GS.f =  (G.f - GS.f).^2
    #print((maximum(GS.f), minimum(GS.f)))
    histogram(GS.f, leg = false, nbins = 100)
    current()
end

###
frac = 1/4
p = Int64(ceil(frac*G.length))

### original
GraphSig_Plot(G, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0., clim = color_limit)
Plots.savefig("original.pdf")

### haar
top_vectors_synthesis(p, dvec_haar, BS_haar, GP, G)
Plots.savefig("syn_haar.pdf")

top_vectors_residual(p, dvec_haar, BS_haar, GP, G)
Plots.savefig("syn_haar_residual_squared.pdf")

#top_vectors_residual_hist(p, dvec_haar, BS_haar, GP, G)
#Plots.savefig("syn_haar_residual_squared_hist.pdf")


### walsh
top_vectors_synthesis(p, dvec_walsh, BS_walsh, GP, G)
Plots.savefig("syn_walsh.pdf")

top_vectors_residual(p, dvec_walsh, BS_walsh, GP, G)
Plots.savefig("syn_walsh_residual_squared.pdf")

#top_vectors_residual_hist(p, dvec_walsh, BS_walsh, GP, G)
#Plots.savefig("syn_walsh_residual_squared_hist.pdf")


### c2f
top_vectors_synthesis(p, dvec_c2f, BS_c2f, GP, G)
Plots.savefig("syn_c2f.pdf")

top_vectors_residual(p, dvec_c2f, BS_c2f, GP, G)
Plots.savefig("syn_c2f_residual_squared.pdf")

#top_vectors_residual_hist(p, dvec_c2f, BS_c2f, GP, G)
#Plots.savefig("syn_c2f_residual_squared_hist.pdf")


### f2c
top_vectors_synthesis(p, dvec_f2c, BS_f2c, GP, G)
Plots.savefig("syn_f2c.pdf")

top_vectors_residual(p, dvec_f2c, BS_f2c, GP, G)
Plots.savefig("syn_f2c_residual_squared.pdf")

#top_vectors_residual_hist(p, dvec_f2c, BS_f2c, GP, G)
#Plots.savefig("syn_f2c_residual_squared_hist.pdf")


### eghwt
top_vectors_synthesis(p, dvec_eghwt, BS_eghwt, GP, G)
Plots.savefig("syn_eghwt.pdf")

top_vectors_residual(p, dvec_eghwt, BS_eghwt, GP, G)
Plots.savefig("syn_eghwt_residual_squared.pdf")

#top_vectors_residual_hist(p, dvec_eghwt, BS_eghwt, GP, G)
#Plots.savefig("syn_eghwt_residual_squared_hist.pdf")



########################################################################
############## Visuallization of Top basis vector#######################
########################################################################

### function to plot the first i-th vector
### Note that here only one vector is plotted instead of multiple vectors.
function top_vectors_plot(dvec::Array{Float64, 2}, BS::BasisSpec, GP::GraphPart, G::GraphSig, i::Int64)
    sorted_ind = sortperm(abs.(dvec[:]), rev = true);
    dvecT = fill(0., size(dvec))
    #dvecT[sorted_ind[i]] = dvec_ghwt[sorted_ind[i]]
    dvecT[sorted_ind[i]] = 1
    f, GS = ghwt_synthesis(dvecT, GP, BS, G)
    GS.name = "The "*string(i)*"-th vector"
    GraphSig_Plot(GS, linewidth = 1., markersize = 4., markercolor = :viridis, markerstrokealpha =0.)
    current()
end
########################## plot top i-th vector
i = 1
### haar
top_vectors_plot(dvec_haar, BS_haar, GP, G, i)
### walsh
top_vectors_plot(dvec_walsh, BS_walsh, GP, G, i)
### c2f
top_vectors_plot(dvec_c2f, BS_c2f, GP, G, i)
### f2c
top_vectors_plot(dvec_f2c, BS_f2c, GP, G, i)
### eghwt
top_vectors_plot(dvec_eghwt, BS_eghwt, GP, G, i)
