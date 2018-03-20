using MTSG, MAT

tmp = matread("MN_MutGauss.mat")
tmp1 = tmp["G"]
G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"])

G = Adj2InvEuc(G)
GP = partition_tree_fiedler(G,:Lrw)
GN = AddNoise(G, SNR = 5.0, noisetype = "gaussian")
dmatrix = ghwt_analysis!(GN, GP=GP)

# through the old way
dvec,BS = ghwt_bestbasis(dmatrix, GP,cfspec=1)
(dvecT,kept) = dvec_Threshold(dvec,"s",0.11,GP,BS)
(f, GS) = ghwt_synthesis(reshape(dvecT,(size(dvecT)[1],1)), GP, BS, GN)
println("The GHWT reconstructed signal has SNR ", snr(G,GS)," dB")

# through time-frequency analysis
bestbasis, bestbasis_tag = ghwt_tf_bestbasis(dmatrix[:,:,1], GP)
bestbasis_T = tf_threshold(bestbasis, GP, 0.11, "s")
(f_tf, GS_tf) = tf_synthesis(bestbasis_T, bestbasis_tag, GP, GN)
println("The tf-analysis GHWT reconstructed signal has SNR ", snr(G,GS_tf)," dB")
