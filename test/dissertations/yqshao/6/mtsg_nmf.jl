#Table 6.2
#Run uility functions at bottom
#Change the setting of Lambda if necessary
using NMF, LinearAlgebra, eghwt, ScikitLearn
@sk_import linear_model:Lasso

k = 5
m = 125
n = 25
sigma = 0.5 #noise parameter
#methodid = 1 corresponds to ALSPGRAD.
#methodid = 2 corresponds to HALS.
methodid = 2


num_pos_eghwt_svd = 0
num_neg_eghwt_svd = 0
num_pos_eghwt_rand = 0
num_neg_eghwt_rand = 0
num_pos_svd_rand = 0
num_neg_svd_rand = 0
num_pos_svd_newsvd = 0
num_neg_svd_newsvd = 0
num_pos_eghwt_newsvd = 0
num_neg_eghwt_newsvd = 0
for i = 1:50
    Wtrue = rand(m,k)
    Htrue = rand(k,n)
    X = Wtrue*Htrue + sigma*rand(m,n)

    #X = rand(m,n)

    ################
    ######################################

    Wtilda, Htilda = ghwtinit(X, k)
    if Wtilda == false
        continue
    end
    Wsvd, Hsvd = NMF.nndsvd(X, k)
    Wrand, Hrand = NMF.randinit(X, k)
    Wnewsvd, Hnewsvd = improved_NNDSVD(X, k)

    iters = 100000
    method = [NMF.ALSPGrad{Float64}(maxiter=iters, tolg=1.0e-6), NMF.CoordinateDescent{Float64}(maxiter=iters, α=0.5, l₁ratio=0.5)]

    a = NMF.solve!(method[methodid], X, copy(Wtilda), copy(Htilda))
    b = NMF.solve!(method[methodid], X, copy(Wsvd), copy(Hsvd))
    c = NMF.solve!(method[methodid], X, copy(Wrand), copy(Hrand))
    d = NMF.solve!(method[methodid], X, copy(Wnewsvd), copy(Hnewsvd))

    #println("initial norm ",norm(X - Wtilda*Htilda)," ",norm(X-Wsvd*Hsvd)," ",norm(X - Wrand*Hrand)," ", norm(X - Wnewsvd*Hnewsvd))
    #println(" ",a.niters," ",b.niters," ",c.niters," ",d.niters)
    #println("final norm",norm(X - a.W*a.H)," ",norm(X-b.W*b.H)," ",norm(X - c.W*c.H)," ", norm(X - d.W*d.H))
    if a.niters < b.niters
        global num_pos_eghwt_svd += 1
    else
        global num_neg_eghwt_svd += 1
    end
    if a.niters < c.niters
        global num_pos_eghwt_rand += 1
    else
        global num_neg_eghwt_rand += 1
    end
    if b.niters < c.niters
        global num_pos_svd_rand += 1
    else
        global num_neg_svd_rand += 1
    end
    if a.niters < d.niters
        global num_pos_eghwt_newsvd += 1
    else
        global num_neg_eghwt_newsvd += 1
    end
    if b.niters < d.niters
        global num_pos_svd_newsvd += 1
    else
        global num_neg_svd_newsvd += 1
    end
end

println("compare eGHWT with NNDSVD ", num_pos_eghwt_svd/(num_pos_eghwt_svd + num_neg_eghwt_svd))
println("compare eGHWT with Random ", num_pos_eghwt_rand/(num_pos_eghwt_rand + num_neg_eghwt_rand))
println("compare eGHWT with NNSVD-LRC ", num_pos_eghwt_newsvd/(num_pos_eghwt_newsvd + num_neg_eghwt_newsvd))




###########################
########################### utility functions
function create_scaling_vectors(GProws::GraphPart)
    GP = GProws
    ### creating scaling vectors
    temp = GP.tag.== 0;
    scaling_index = findall(temp);# index of the scaling vectors .i.e. vectors with tag == 0

    num_scaling = length(scaling_index);
    scaling_vectors = fill(0., (size(GP.tag,1), num_scaling));

    for i = 1:num_scaling
        ind = scaling_index[i]
        dvec = fill(0., (size(GP.tag,1),1))
        dvec[ind[1],1] = 1.;
        BS = bs_level(GP, ind[2] - 1)
        scaling_vectors[:,i] = ghwt_synthesis(dvec, GP, BS)
    end
    return scaling_vectors
end

function create_PHI(scaling_vectors_rows, scaling_vectors_cols, m, n)
    PHI_right = repeat(scaling_vectors_cols[:,1:n], outer = (1, m))
    PHI_left = repeat(scaling_vectors_rows[:,1:m], inner = (1, n))
    return PHI_left, PHI_right
end

function ghwtinit(matrix, k)
    GProws, GPcols = partition_tree_matrixDhillon(matrix)
    dmatrix = ghwt_analysis_2d(matrix, GProws, GPcols)
    scaling_vectors_rows = create_scaling_vectors(GProws)
    scaling_vectors_cols = create_scaling_vectors(GPcols)

    PHI_left, PHI_right = create_PHI(scaling_vectors_rows, scaling_vectors_cols, 10, 10)
    scaling_vectors = fill(0., (length(matrix[:]), size(PHI_left,2)))
    for i = 1:size(PHI_left,2)
        scaling_vectors[:,i] = (PHI_left[:,i]*PHI_right[:,i]')[:]
    end
    ##########################
    ####################################### perform lasso

    f = matrix[:]
    #Lambda = range(0.00001,stop = 0.0001, length = 10)
    Lambda = range(0.001,stop = 0.01, length = 10)# m=125,n=25,k=5,sigma = 0.5,
    #Lambda = range(0.001, stop = 0.01, length = 10)
    B = fill(0., (size(scaling_vectors,2), length(Lambda))) ### initialize coefficients
    for i = 1:length(Lambda)
        lasso = Lasso(alpha = Lambda[i], positive = true, fit_intercept = false)
        fit!(lasso, scaling_vectors, f)
        B[:,i] = lasso.coef_
    end

    j = 1
    if sum((B[:,j].!=0)) < k
        return false, false
    else
        while j <10 && sum(B[:,j+1].!=0) >= k
            j += 1
        end
    end
    index = findall(B[:,j].!=0)[1:k]
    rescaleW = repeat(sqrt.(B[index,j]'),size(PHI_left,1),1)
    rescaleH = repeat(sqrt.(B[index,j]'),size(PHI_right,1),1)
    Wtilda = PHI_left[:,index].*rescaleW
    Htilda = Array{Float64,2}((PHI_right[:,index].*rescaleH)')
    return Wtilda, Htilda
end


function improved_NNDSVD(matrix, k)
    m,n = size(matrix)
    U,S,V = svd(matrix)
    W = zeros(m,k)
    H = zeros(k,n)
    Y = U*sqrt.(Diagonal(S))
    Z = sqrt.(Diagonal(S))*transpose(V)

    W[:,1] = abs.(Y[:,1])
    H[1,:] = abs.(Z[1,:])
    j = 2
    for i = 2:k
        if mod(i,2) == 0
            W[:,i] = Y[:,j]
            H[i,:] = Z[i,:]
        else
            W[:,i] = .-Y[:,j]
            H[i,:] = .-Z[i,:]
            j += 1
        end
        W[W[:,i] .< 0,i] .= 0
        H[i,H[i,:] .< 0] .= 0
    end
    return W, H
end
