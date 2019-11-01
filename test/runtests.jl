# This is a very preliminary test function;
# Just small scale examples, e.g., P6 with 10 random signals. More coming!
using MTSG, LinearAlgebra, SparseArrays, JLD2
##############################################################
# 1. Testing basic GHWT functions on 10 random signals on P6 #
##############################################################
println("1. Testing basic GHWT functions")
JLD2.@load "runtests_data/path6randn10.jld2" G
G = gpath(6, G["tmp"])
GP = partition_tree_fiedler(G)
dc2f = ghwt_analysis!(G, GP=GP)
# Check the Haar basis
BH=bs_haar(GP)
dvec=dmatrix2dvec(dc2f, GP, BH)
frecon=ghwt_synthesis(dvec,GP,BH)
println("Relative L2 error of the Haar transform: ", norm(G.f-frecon)/norm(G.f))

# Check the Best Basis
bbc2f = ghwt_c2f_bestbasis(dc2f, GP)
levlist = collect(enumerate([4; 4; 3; 3; 3; 3]))
println("The true BB levlist: ", levlist)
println("The comp BB levlist: ", bbc2f[2].levlist)
println("They are equal: ", levlist == bbc2f[2].levlist)
# levlengths = collect(enumerate([1; 1; 1; 2; 1]))
# println("The true BB levlengths: ", levlengths)
# println("The comp BB levlengths: ", bbc2f[2].levlengths)
JLD2.@load "runtests_data/bbcoef.jld2" tmp2
println("The relative L2 error of the BB coefs: ", norm(tmp2["bbcoef"]-bbc2f[1])/norm(tmp2["bbcoef"]))
println("\n")


###############################################################################
# 2. Testing time-frequency adapted GHWT functions on smoothing the Minnesota #
#   roadmap signal (Temporarily skipped)                                      #
###############################################################################
# println("2. Testing time-frequency adapted GHWT functions on smoothing the minnesota roadmap signal")
# JLD2.@load "runtests_data/MN_MutGauss.jld2" G
# tmp1 = G["G"]
# G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])
#
# G = Adj2InvEuc(G)
# GP = partition_tree_fiedler(G,:Lrw)
# GN = AddNoise(G, SNR = 5.0, noisetype = "gaussian")
# dmatrix = ghwt_analysis!(GN, GP=GP)
#
# println("The noisy signal has SNR 5.0")
# # through the old way
# dvec,BS = ghwt_bestbasis(dmatrix, GP,cfspec=1)
# (dvecT,kept) = dvec_Threshold(dvec,"s",0.11,GP,BS)
# (f, GS) = ghwt_synthesis(reshape(dvecT,(size(dvecT)[1],1)), GP, BS, GN)
# println("The GHWT smoothed signal has SNR ", snr(G,GS)," dB")
#
# # through time-frequency analysis
# dvec_tf, BS_tf = ghwt_tf_bestbasis(dmatrix, GP)
# (dvecT_tf,kept) = dvec_Threshold(dvec_tf,"s",0.11,GP,BS_tf)
# (f_tf, GS_tf) = ghwt_synthesis(reshape(dvecT_tf,(size(dvecT_tf)[1],1)), GP, BS_tf, GN)
# println("The tf-analysis GHWT smoothed signal has SNR ", snr(G,GS_tf)," dB")
#############################

#############################################################################
# 2. Testing time-frequency adapted GHWT functions on a simple signal on P6 #
#############################################################################
f = Array{Float64}([2. -2. 1. 3. -1. -2.]')
#G = GraphSig(SparseMatrixCSC(diagm(1 => ones(5))),f = f)
G = gpath(6, f)
GP = partition_tree_fiedler(G,:Lrw)
dmatrix = ghwt_analysis!(G, GP=GP)
println("2. Testing time-frequency adapted GHWT functions on path signal: ", f)
println("The original signal has L1 norm: ", norm(f,1))

# through the old way
dvec,BS = ghwt_bestbasis(dmatrix, GP, cfspec=1.)
(f, GS) = ghwt_synthesis(dvec, GP, BS, G)
println("The coefficient vectors of GHWT best basis has L1 norm: ", norm(dvec,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f-f)/norm(G.f))

# through the old way but restricted to level j_start to level j_end in c2f dictionary
j_start = 3
j_end = 4
dvec_restricted,BS_restricted = ghwt_bestbasis(dmatrix, GP, cfspec=1., j_start = j_start, j_end = j_end)
(f_restricted, GS_restricted) = ghwt_synthesis(dvec_restricted, GP, BS_restricted, G)
println("The coefficient vectors of GHWT best basis restricted from level $(j_start) to level $(j_end) in c2f dictionary has L1 norm: ", norm(dvec_restricted,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f-f_restricted)/norm(G.f))

# through the time-frequency adapted analysis
dvec_tf, BS_tf = ghwt_tf_bestbasis(dmatrix, GP, cfspec=1.)
(f_tf, GS_tf) = ghwt_synthesis(dvec_tf, GP, BS_tf, G)
println("The coefficient vectors of time-frequency adapted GHWT best basis has L1 norm: ", norm(dvec_tf,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f-f_tf)/norm(G.f))
println("\n")


#################################################################################
# 3. Testing time-frequency adapted 2D GHWT functions using a simple 3x3 matrix #
#################################################################################
println("3. Testing time-frequency adapted 2D GHWT functions")
matrix = [ 1.0 2.0 3.0; 4.0 5.0 6.0]
#matrix = matread("termdoc.mat")["matrix"]

# expand matrix in 2 directions (rows and cols)
GProws, GPcols = partition_tree_matrixDhillon(matrix)
dmatrix = ghwt_analysis_2d(matrix, GProws, GPcols)

# find the best basis using the time-frequency analysis
# infovec indicate the location of each coefficient in dmatrix
dvec_tf, infovec_tf = ghwt_tf_bestbasis_2d(matrix, GProws, GPcols)
matrix_tf = ghwt_tf_synthesis_2d(dvec_tf, infovec_tf, GProws, GPcols)

# find the best basis using the GHWT
dvec_ghwt, BSrows, BScols = ghwt_bestbasis_2d(matrix, GProws, GPcols)
matrix_ghwt = ghwt_synthesis_2d(dvec_ghwt, GProws, GPcols, BSrows, BScols)

println("The toy matrix [1 2 3; 4 5 6] has 1-vecnorm as ", norm(matrix,1))
println("The eghwt bestbasis of toy matrix [1 2 3; 4 5 6] has 1-vecnorm as ", norm(dvec_tf,1))
println("The ghwt bestbasis of toy matrix [1 2 3; 4 5 6] has 1-vecnorm as ", norm(dvec_ghwt,1))
println("Relative L2 error of the eghwt synthesized matrix: ", norm(matrix - matrix_tf, 2)/norm(matrix,2))
println("Relative L2 error of the ghwt synthesized matrix: ", norm(matrix - matrix_ghwt, 2)/norm(matrix,2))
println("\n")

#############################################################################
# 4. Testing HGLET functions and hybrid methods related functions on RGC100 #
#############################################################################
println("4. Testing HGLET functions and hybrid methods related functions")
#JLD2.@load "runtests_data/Dendrite.jld2" G
#G = G["G"]
# convert to SparseMatrixCSC{Float64,Int} where Int is machine dependent.
#A = Array(G["W"])
#B = sparse(A)
#G = GraphSig(B, xy = G["xy"], f = G["f"], name = G["name"])

JLD2.@load "runtests_data/Toronto_new.jld2"
tmp1 = toronto["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])


GP = partition_tree_fiedler(G)

dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G, GP) # expansion coefficients of 3-way HGLET bases
dmatrixG = ghwt_analysis!(G,GP = GP) # expansion coefficients of GHWT

dvec5,BS5,trans5,p5 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH, dmatrixG = dmatrixG,
dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym) # best-basis among all combinations of bases

fS5, GS5 = HGLET_GHWT_Synthesis(reshape(dvec5,(size(dvec5)[1],1)),GP,BS5,trans5,G)

println("The original signal has L1 norm: ", norm(G.f,1))
println("The coefficients of best-basis selected from hybrid method has L1 norm: ", norm(dvec5,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f-fS5)/norm(G.f))
println("\n")

###########################################################################
# 5. Testing GraphSig_Plot on mutilated Gaussian signal on the MN roadmap #
###########################################################################
println("5. GraphSig_Plot on mutilated Gaussian signal on the MN roadmap")
println("... this may not display the plot if you run it via the package test.")
JLD2.@load "runtests_data/MN_MutGauss.jld2" G
tmp1 = G["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])
G = Adj2InvEuc(G)
GraphSig_Plot(G)

# End of runtests.jl
