#using MTSG, MAT

## 0.Preliminaries
#include("src\\MTSG.jl")
using MTSG, MAT
#G = matread("test\\Dendrite.mat")
G = matread("Dendrite.mat")

# paratition the graph
G = G["G"]
G = GraphSig(G["W"], xy = G["xy"], f = G["f"], name = G["name"])
GP = partition_tree_fiedler(G)


## 1. Analyze the signal
dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G, GP) # expansion coefficients of 3-way HGLET bases
dmatrixG = ghwt_analysis!(G,GP = GP) # expansion coefficients of GHWT

# 1) HGLET best basis with L
dvec1,BS1,_,p1 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH)
r1 = orth2relerror(dvec1)

# 2) Laplacian eigenvectors
BS2 = LevelBasisSpec(GP,0)
dvec2 = dmatrix2dvec(dmatrixH,GP,BS2)
r2 = orth2relerror(dvec2[:,1])

# 3) GHWT best basis
dvec3,BS3,_,p3 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixG = dmatrixG)
r3 = orth2relerror(dvec3)

# 4) Haar basis
BS4 = bs_haar(GP)
dvec4 = dmatrix2dvec(dmatrixG, GP, BS4)
r4 = orth2relerror(dvec4[:,1])

# 5) Hybrid best basis (HGLET with L/Lrw/Lsym + GHWT)
dvec5,BS5,trans5,p5 = HGLET_GHWT_BestBasis_minrelerror(GP,G,dmatrixH = dmatrixH, dmatrixG = dmatrixG,
dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym) # best-basis among all combinations of bases
B5,_ = HGLET_GHWT_Synthesis(eye(N),GP,BS5,trans5,G)
r5 = nonorth2relerror(dvec5[:,1],B5)


# 6) Walsh basis
BS6 = LevelBasisSpec(GP,0)
dvec6 = dmatrix2dvec(dmatrixG,GP,BS6)
r6 = orth2relerror(dvec6[:,1])
