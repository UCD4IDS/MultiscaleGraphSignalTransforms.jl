using MTSG, MAT

G = matread("Dendrite.mat")
G = G["G"]
G = GraphSig(G["W"], xy = G["xy"], f = G["f"], name = G["name"])
GP = partition_tree_fiedler(G)

dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G, GP) # expansion coefficients of 3-way HGLET bases
dmatrixG = ghwt_analysis!(G,GP = GP) # expansion coefficients of GHWT

dvec5,BS5,trans5,p5 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH, dmatrixG = dmatrixG,
dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym) # best-basis among all combinations of bases
