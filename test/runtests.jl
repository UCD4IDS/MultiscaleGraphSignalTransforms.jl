# This is a very preliminary test function; just a copy of bbtest.jl of small scale P6 with 10 random signals. More coming!
using MTSG, MAT
###########################################################################################
# Testing basic GHWT functions #
###########################################################################################

println("Testing basic GHWT functions")
tmp=matread("path6randn10.mat")
G=gpath(6, tmp["tmp"])
GP=partition_tree_fiedler(G)
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
tmp2=matread("bbcoef.mat")
println("The relative L2 error of the BB coefs: ", norm(tmp2["bbcoef"][:]-bbc2f[1][:])/norm(tmp2["bbcoef"]))


###########################################################################################
# Testing time-frequency adapted GHWT functions on smoothing the minnesota roadmap signal #
###########################################################################################
println("Testing time-frequency adapted GHWT functions on smoothing the minnesota roadmap signal")
tmp = matread("MN_MutGauss.mat")
tmp1 = tmp["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])

G = Adj2InvEuc(G)
GP = partition_tree_fiedler(G,:Lrw)
GN = AddNoise(G, SNR = 5.0, noisetype = "gaussian")
dmatrix = ghwt_analysis!(GN, GP=GP)

println("The noisy signal has SNR 5.0")
# through the old way
dvec,BS = ghwt_bestbasis(dmatrix, GP,cfspec=1)
(dvecT,kept) = dvec_Threshold(dvec,"s",0.11,GP,BS)
(f, GS) = ghwt_synthesis(reshape(dvecT,(size(dvecT)[1],1)), GP, BS, GN)
println("The GHWT smoothed signal has SNR ", snr(G,GS)," dB")

# through time-frequency analysis
bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix[:,:,1], GP)
bestbasis_T = tf_threshold(bestbasis, GP, 0.11, "s")
(f_tf, GS_tf) = tf_synthesis(bestbasis_T, bestbasis_tag, GP, GN)
println("The tf-analysis GHWT smoothed signal has SNR ", snr(G,GS_tf)," dB")


###########################################################################################
# Testing time-frequency adapted 2d GHWT functions #
###########################################################################################
println("Testing time-frequency adapted 2d GHWT functions")
matrix = [ 1.0 2.0 3.0; 4.0 5.0 6.0]
#matrix = matread("termdoc.mat")["matrix"]


# expand matrix in 2 directions (rows and cols)
dmatrix, GProws, GPcols = ghwt_tf_init_2d(matrix)

# find the best basis using the time-frequency analysis
# infovec indicate the location of each coefficient in dmatrix
Bbasis, infovec = ghwt_tf_bestbasis_2d(dmatrix, GProws, GPcols)

println("The bestbasis of toy matrix [1, 2, 3; 4, 5, 6] has Frobenius norm as ", vecnorm(Bbasis,1))

###########################################################################################
# Testing HGLET functions and hybrid methods related functions #
###########################################################################################
println("Testing HGLET functions and hybrid methods related functions")
G = matread("Dendrite.mat")
G = G["G"]
G = GraphSig(G["W"], xy = G["xy"], f = G["f"], name = G["name"])
GP = partition_tree_fiedler(G)

dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G, GP) # expansion coefficients of 3-way HGLET bases
dmatrixG = ghwt_analysis!(G,GP = GP) # expansion coefficients of GHWT

dvec5,BS5,trans5,p5 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH, dmatrixG = dmatrixG,
dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym) # best-basis among all combinations of bases

println("The bestbasis of hybrid methods of MN data has Frobenius norm as ", vecnorm(dvec5,1))
