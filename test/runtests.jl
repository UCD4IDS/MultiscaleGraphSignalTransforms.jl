include("..\\src\\MTSG.jl")
# This is a very preliminary test function; just a copy of bbtest.jl of small scale P6 with 10 random signals. More coming!
using MTSG, MAT
tmp=matread("test\\path6randn10.mat")
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
levlist = Vector{UInt8}([4, 4, 3, 3, 3])
println("The true BB levlist: ", levlist')
println("The comp BB levlist: ", (bbc2f[2].levlist)')
levlengths = Vector{UInt8}([1, 1, 1, 2, 1])
println("The true BB levlengths: ", levlengths')
println("The comp BB levlengths: ", (bbc2f[2].levlengths)')
tmp2=matread("test\\bbcoef.mat")
println("The relative L2 error of the BB coefs: ", norm(tmp2["bbcoef"][:]-bbc2f[1][:])/norm(tmp2["bbcoef"]))
