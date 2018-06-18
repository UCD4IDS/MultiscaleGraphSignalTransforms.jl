using MTSG, MAT


include("..\\src\\GHWT_tf_2d.jl")


matrix = [ 1.0 2.0 3.0; 4.0 5.0 6.0]
#matrix = matread("termdoc.mat")["matrix"]


# expand matrix in 2 directions (rows and cols)
dmatrix, GProws, GPcols = ghwt_tf_init_2d(matrix)

# find the best basis using the time-frequency analysis
# infovec indicate the location of each coefficient in dmatrix
Bbasis, infovec = ghwt_tf_bestbasis_2d(dmatrix, GProws, GPcols)
