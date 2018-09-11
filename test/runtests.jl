# This is a very preliminary test function; just a copy of bbtest.jl of small scale P6 with 10 random signals. More coming!
Pkg.add("JLD2")
using MTSG, LinearAlgebra, SparseArrays, JLD2
###########################################################################################
# Testing basic GHWT functions #
###########################################################################################

println("1. Testing basic GHWT functions")
@load "path6randn10.jld2" G
G = gpath(6, G["tmp"])
GP = partition_tree_fiedler(G)
dc2f = ghwt_analysis!(G, GP=GP)
# Check the Haar basis
BH=bs_haar(GP)
dvec=dmatrix2dvec(dc2f, GP, BH)
frecon=ghwt_synthesis(dvec,GP,BH)
println("Relative L2 error of the Haar transform: ", norm(G.f[:]-frecon[:])/norm(G.f[:]))

# Check the Best Basis
bbc2f = ghwt_c2f_bestbasis(dc2f, GP)
levlist = Vector{Int}([4, 4, 3, 3, 3])
println("The true BB levlist: ", levlist')
println("The comp BB levlist: ", (bbc2f[2].levlist)')
levlengths = Vector{Int}([1, 1, 1, 2, 1])
println("The true BB levlengths: ", levlengths')
println("The comp BB levlengths: ", (bbc2f[2].levlengths)')
@load "bbcoef.jld2" tmp2
println("The relative L2 error of the BB coefs: ", norm(tmp2["bbcoef"][:]-bbc2f[1][:])/norm(tmp2["bbcoef"]))
println("\n")


###########################################################################################
# Testing time-frequency adapted GHWT functions on smoothing the minnesota roadmap signal #
###########################################################################################
# println("2. Testing time-frequency adapted GHWT functions on smoothing the minnesota roadmap signal")
# tmp = matread("MN_MutGauss.mat")
# tmp1 = tmp["G"]
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
# bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix[:,:,1], GP)
# bestbasis_T = tf_threshold(bestbasis, GP, 0.11, "s")
# (f_tf, GS_tf) = tf_synthesis(bestbasis_T, bestbasis_tag, GP, GN)
# println("The tf-analysis GHWT smoothed signal has SNR ", snr(G,GS_tf)," dB")


#############################
f = Array{Float64}([2. -2. 1. 3. -1. -2.]')
G = GraphSig(SparseMatrixCSC(diagm(1 => ones(5))),f = f)
GP = partition_tree_fiedler(G,:Lrw)
dmatrix = ghwt_analysis!(G, GP=GP)
println("2. Testing time-frequency adapted GHWT functions on path signal: ", f)
println("The original signal has L1 norm: ", norm(f,1))

# through the old way
dvec,BS = ghwt_bestbasis(dmatrix, GP,cfspec=1)
(f, GS) = ghwt_synthesis(reshape(dvec,(size(dvec)[1],1)), GP, BS, G)
println("The coefficient vectors of GHWT best basis has L1 norm: ", norm(dvec,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f[:]-f[:])/norm(G.f[:]))


# through the time-frequency analysis
bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix[:,:,1], GP)
(f_tf, GS_tf) = tf_synthesis(bestbasis, bestbasis_tag, GP, G)
println("The coefficient vectors of time-frequency adapted GHWT best basis has L1 norm: ", norm(bestbasis,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f[:]-f_tf[:])/norm(G.f[:]))

println("\n")








###########################################################################################
# Testing time-frequency adapted 2d GHWT functions #
###########################################################################################
println("3. Testing time-frequency adapted 2d GHWT functions")
matrix = [ 1.0 2.0 3.0; 4.0 5.0 6.0]
#matrix = matread("termdoc.mat")["matrix"]


# expand matrix in 2 directions (rows and cols)
dmatrix, GProws, GPcols = ghwt_tf_init_2d(matrix)
dmatrix = Array{Float64,2}(dmatrix)

# find the best basis using the time-frequency analysis
# infovec indicate the location of each coefficient in dmatrix
Bbasis, infovec = ghwt_tf_bestbasis_2d(dmatrix, GProws, GPcols)

BBmatrix = zeros(size(dmatrix))
for i in 1:size(infovec,2)
    BBmatrix[infovec[1,i],infovec[2,i]] = Bbasis[i]
end

matrix_r_temp = ghwt_synthesis_2d(BBmatrix, GProws, GPcols)
matrix_r = zeros(size(matrix))
matrix_r[GProws.ind,GPcols.ind] = matrix_r_temp

println("The toy matrix [1, 2, 3; 4, 5, 6] has 1-vecnorm as ", norm(matrix,1))
println("The bestbasis of toy matrix [1, 2, 3; 4, 5, 6] has 1-vecnorm as ", norm(Bbasis,1))
println("Relative L2 error of the synthesized matrix: ", norm(matrix[:] - matrix_r[:], 2)/norm(matrix[:],2))

println("\n")

###########################################################################################
# Testing HGLET functions and hybrid methods related functions #
###########################################################################################
println("4. Testing HGLET functions and hybrid methods related functions")
@load "Dendrite.jld2" G
G = G["G"]
G = GraphSig(G["W"], xy = G["xy"], f = G["f"], name = G["name"])
GP = partition_tree_fiedler(G)

dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G, GP) # expansion coefficients of 3-way HGLET bases
dmatrixG = ghwt_analysis!(G,GP = GP) # expansion coefficients of GHWT

dvec5,BS5,trans5,p5 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH, dmatrixG = dmatrixG,
dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym) # best-basis among all combinations of bases

fS5, GS5 = HGLET_GHWT_Synthesis(reshape(dvec5,(size(dvec5)[1],1)),GP,BS5,trans5,G)

println("The original signal has L1 norm: ", norm(G.f,1))
println("The coefficients of best-basis selected from hybrid method has L1 norm: ", norm(dvec5,1))
println("Relative L2 error of the synthesized signal: ", norm(G.f[:]-fS5[:])/norm(G.f[:]))

println("\n")
