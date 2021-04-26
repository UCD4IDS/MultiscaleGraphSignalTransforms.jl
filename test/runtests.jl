# This is a very preliminary test function;
# Just small scale examples, e.g., P6 with 10 random signals. More coming!
using Test, MultiscaleGraphSignalTransforms, LinearAlgebra, SparseArrays, JLD2, Plots
@testset "MultiscalGraphSignalTransforms.jl" begin
##############################################################
# 1. Testing basic GHWT functions on 10 random signals on P6 #
##############################################################
    @testset "1. Testing basic GHWT functions" begin
        println("1. Testing basic GHWT functions")
        JLD2.@load "runtests_data/path6randn10.jld2" G
        G = gpath(6, G["tmp"])
        GP = partition_tree_fiedler(G)
        dc2f = ghwt_analysis!(G, GP=GP)

        @testset "Check the Haar basis" begin
            BH=bs_haar(GP)
            dvec=dmatrix2dvec(dc2f, GP, BH)
            frecon=ghwt_synthesis(dvec,GP,BH)
            println("Relative L2 error of the Haar transform: ", norm(G.f-frecon)/norm(G.f))
            @test norm(G.f-frecon)/norm(G.f) < 10 * eps()
        end

        @testset "Check the best basis and its expansion coefficients" begin
            bbc2f = ghwt_c2f_bestbasis(dc2f, GP)
            levlist = collect(enumerate([4; 4; 3; 3; 3; 3]))
            println("The true BB levlist: ", levlist)
            println("The comp BB levlist: ", bbc2f[2].levlist)
            println("They are equal: ", levlist == bbc2f[2].levlist)
            @test levlist == bbc2f[2].levlist
            JLD2.@load "runtests_data/bbcoef.jld2" tmp2
            println("The relative L2 error of the BB coefs: ", norm(tmp2["bbcoef"]-bbc2f[1])/norm(tmp2["bbcoef"]))
            @test norm(tmp2["bbcoef"]-bbc2f[1])/norm(tmp2["bbcoef"]) < 10 * eps()
        end
        println("\n")
    end

################################################################
# 2. Testing GHWT/eGHWT functions on a simple synthetic signal on P6 #
################################################################
    @testset "2. GHWT/eGHWT tests on P6" begin
        f = Array{Float64}([2. -2. 1. 3. -1. -2.]')
        G = gpath(6, f)
        GP = partition_tree_fiedler(G; method = :Lrw)
        dmatrix = ghwt_analysis!(G, GP=GP)
        println("2. Testing the GHWT functions on a path signal: ", f)
        println("The original signal has L1 norm: ", norm(f,1))

        @testset "The full level GHWT best basis" begin
            dvec,BS = ghwt_bestbasis(dmatrix, GP, cfspec=1.)
            (f, GS) = ghwt_synthesis(dvec, GP, BS, G)
            println("The coefficient vectors of the GHWT best basis has L1 norm: ", norm(dvec,1))
            println("Relative L2 error of the synthesized signal: ", norm(G.f-f)/norm(G.f))
            @test norm(G.f-f)/norm(G.f) < 10 * eps()
        end

        @testset "The GHWT best basis restricted levels" begin
            j_start = 3
            j_end = 4
            dvec_restricted,BS_restricted = ghwt_bestbasis(dmatrix, GP, cfspec=1., j_start = j_start, j_end = j_end)
            (f_restricted, GS_restricted) = ghwt_synthesis(dvec_restricted, GP, BS_restricted, G)
            println("The coefficient vectors of the GHWT best basis restricted from level $(j_start) to level $(j_end) in c2f dictionary has L1 norm: ", norm(dvec_restricted,1))
            println("Relative L2 error of the synthesized signal: ", norm(G.f-f_restricted)/norm(G.f))
            @test norm(G.f-f_restricted)/norm(G.f) < 10 * eps()
        end

        @testset "The eGHWT best basis" begin
            dvec_tf, BS_tf = ghwt_tf_bestbasis(dmatrix, GP, cfspec=1.)
            (f_tf, GS_tf) = ghwt_synthesis(dvec_tf, GP, BS_tf, G)
            println("The coefficient vectors of the time-frequency adapted GHWT best basis has L1 norm: ", norm(dvec_tf,1))
            println("Relative L2 error of the synthesized signal: ", norm(G.f-f_tf)/norm(G.f))
            @test norm(G.f-f_tf)/norm(G.f) < 10 * eps()
        end
        println("\n")
    end

#################################################################################
# 3. Testing time-frequency adapted 2D GHWT functions using a simple 3x3 matrix #
#################################################################################
    @testset "3. 2D GHWT/eGHWT on 3x3 matrix" begin
        println("3. Testing 2D GHWT/eGHWT functions")
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
        println("The eGHWT best basis of toy matrix [1 2 3; 4 5 6] has 1-vecnorm as ", norm(dvec_tf,1))
        println("The GHWT best basis of toy matrix [1 2 3; 4 5 6] has 1-vecnorm as ", norm(dvec_ghwt,1))
        println("Relative L2 error of the eGHWT synthesized matrix: ", norm(matrix - matrix_tf, 2)/norm(matrix,2))
        @test norm(matrix - matrix_tf, 2)/norm(matrix,2) < 10 * eps()
        println("Relative L2 error of the GHWT synthesized matrix: ", norm(matrix - matrix_ghwt, 2)/norm(matrix,2))
        @test norm(matrix - matrix_ghwt, 2)/norm(matrix,2) < 10 * eps()
        println("\n")
    end

###################################################################
# 4. Testing HGLET functions and hybrid methods on a real dataset #
###################################################################

    @testset "4. HGLET on Toronto Street Network" begin
        println("4. HGLET on Toronto Street Network")
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
        @test norm(G.f-fS5)/norm(G.f) < 10 * eps()
        println("\n")
    end

###########################################################################
# 5. Testing GraphSig_Plot on mutilated Gaussian signal on the MN roadmap #
###########################################################################
    @testset "5. GraphSig_Plot on the MN roadmap" begin
        println("5. GraphSig_Plot on the MN roadmap")
        println("... this may not display the plot if you run it via the package test.")
        JLD2.@load "runtests_data/MN_MutGauss.jld2" G
        tmp1 = G["G"]
        G=GraphSig(tmp1["W"],xy=tmp1["xy"],f=tmp1["f"],name =tmp1["name"],plotspecs = tmp1["plotspecs"])
        G = Adj2InvEuc(G)
        plt = GraphSig_Plot(G)
        display(plt)
        @test typeof(plt) == Plots.Plot{Plots.GRBackend}
        println("\n")
    end

###################################################
# 6. Testing LP-HGLET functions on a real dataset #
###################################################

    @testset "6. LP-HGLET on RGC100" begin
        println("6. LP-HGLET on RGC100")
        JLD2.@load "runtests_data/Dendrite.jld2" G
        G = G["G_3D"]
        f = G["f"]
        G_Sig = GraphSig(1.0 * G["W"]; xy = G["xy"], f = f)
        G_Sig = Adj2InvEuc(G_Sig)
        GP = partition_tree_fiedler(G_Sig; swapRegion = false)

        dmatrixlH, _ = LPHGLET_Analysis_All(G_Sig, GP; ϵ = 0.3)
        dvec_lphglet, BS_lphglet, _ = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixlH)

        fS, GS = LPHGLET_Synthesis(reshape(dvec_lphglet, (length(dvec_lphglet), 1)), GP, BS_lphglet, G_Sig; ϵ = 0.3)
        println("The original signal has L1 norm: ", norm(f, 1))
        println("The coefficients of LP-HGLET best-basis has L1 norm: ", norm(dvec_lphglet, 1))
        println("Relative L2 error of the synthesized signal: ", norm(f - fS) / norm(f))
        @test norm(f - fS) / norm(f) < 20 * eps()
        println("\n")
    end

#########################################################
# 7. Testing VM-NGWP, PC-NGWP, LP-NGWP functions on P64 #
#########################################################

    @testset "7. testing VM-NGWP/PC-NGWP/LP-NGWP on P64" begin
        println("7. testing VM-NGWP/PC-NGWP/LP-NGWP on P64")
        N = 64
        G = gpath(N)
        # use Chebyshev polynomial T₅(x) (x ∈ [0, 1]) as the path signal
        G.f = reshape([16 * x^5 - 20 * x^3 + 5 * x for x in LinRange(0, 1, N)], (N, 1))
        # compute graph Laplacian eigenvectors
        W = G.W
        L = diagm(sum(W; dims = 1)[:]) - W
        𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1,:])'
        # build Gstar
        Gstar = GraphSig(W)
        GP = partition_tree_fiedler(G; swapRegion = false)
        GP_dual = partition_tree_fiedler(Gstar; swapRegion = false)
        GP_primal = pairclustering(𝚽, GP_dual) # for PC-NGWP
        # construct NGWP dictionaries
        VM_NGWP = vm_ngwp(𝚽, GP_dual)
        PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal)
        LP_NGWP = lp_ngwp(𝚽, W, GP_dual; ϵ = 0.3)
        # NGWP analysis
        dmatrix_VM = ngwp_analysis(G, VM_NGWP)
        dvec_vm_ngwp, BS_vm_ngwp = ngwp_bestbasis(dmatrix_VM, GP_dual)
        dmatrix_PC = ngwp_analysis(G, PC_NGWP)
        dvec_pc_ngwp, BS_pc_ngwp = ngwp_bestbasis(dmatrix_PC, GP_dual)
        dmatrix_LP = ngwp_analysis(G, LP_NGWP)
        dvec_lp_ngwp, BS_lp_ngwp = ngwp_bestbasis(dmatrix_LP, GP_dual)

        println("The original signal has L2 norm: ", norm(G.f))
        println("The coefficients of VM-NGWP best-basis has L2 norm: ", norm(dvec_vm_ngwp))
        @test abs(norm(G.f) - norm(dvec_vm_ngwp)) / norm(G.f) < 10 * eps()
        println("The coefficients of PC-NGWP best-basis has L2 norm: ", norm(dvec_pc_ngwp))
        @test abs(norm(G.f) - norm(dvec_pc_ngwp)) / norm(G.f) < 10 * eps()
        println("The coefficients of LP-NGWP best-basis has L2 norm: ", norm(dvec_lp_ngwp))
        @test abs(norm(G.f) - norm(dvec_lp_ngwp)) / norm(G.f) < 10 * eps()
        println("\n")
    end


end
# End of runtests.jl
