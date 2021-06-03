__precompile__()

module MultiscaleGraphSignalTransforms

include("GraphSignal.jl")
include("GraphPartition.jl")
include("BasisSpecification.jl")
include("common.jl")
include("GHWT.jl")
include("GHWT_2d.jl")
include("GHWT_tf_1d.jl")
include("GHWT_tf_2d.jl")
include("HGLET.jl")
include("GraphSig_Plot.jl")
include("gplot.jl")
include("partition_fiedler.jl")

using Reexport
@reexport using .GraphSignal, .GraphPartition, .BasisSpecification, .GHWT, .GHWT_2d, .GHWT_tf_1d, .GHWT_tf_2d, .HGLET

export dvec2dmatrix, dmatrix2dvec, levlist2levlengths!, bsfull, bs_haar, bs_level, bs_walsh, dvec_Threshold, rs_to_region, GraphSig_Plot, gplot, gplot!, partition_fiedler
export cost_functional, dmatrix_flatten

## export functions of NGWP.jl
using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs, Clustering
using JuMP, Clp, Optim, Statistics, QuadGK, Arpack
import Plots: plot, plot!, scatter, scatter!
import StatsBase: crosscor

include("dualgraph.jl")
include("eigDAG_Distance.jl")
include("eigHAD_Distance.jl")
include("eigROT_Distance.jl")
include("eigsROT_Distance.jl")
include("eigTSD_Distance.jl")
include("helpers.jl")
include("LP-HGLET.jl")
include("LP-NGWP.jl")
include("natural_distances.jl")
include("utils.jl")
include("varimax.jl")
include("NGWF.jl")
include("NGWP.jl")
include("PC-NGWP.jl")
include("SunFlowerGraph.jl")
include("VM-NGWP.jl")

export eigDAG_Distance, eigHAD_Distance, eigHAD_Affinity
export eigROT_Distance, ROT_Distance, eigEMD_Distance
export eigsROT_Distance
export eigTSD_Distance, K_functional
export natural_eigdist
export SunFlowerGraph, dualgraph
export pc_ngwp, pairclustering, mgslp
export vm_ngwp, varimax
export lp_ngwp, rising_cutoff, find_pairinds, pair_inds_shadding, lp_ngwp_analysis
export LPHGLET_Synthesis, LPHGLET_Analysis_All, HGLET_dictionary, LPHGLET_dictionary
export unitary_folding_operator, keep_folding!
export ngwp_analysis, ngwp_bestbasis, NGWP_jkl
export natural_eigdist
export nat_spec_filter, ngwf_all_vectors, rngwf_all_vectors, ngwf_vector, frame_approx, rngwf_lx
export scatter_gplot, scatter_gplot!
export standardize_eigenvectors!, spike, characteristic, χ, sort_wavelets, transform2D
export getall_expansioncoeffs, approx_error_plot, getall_expansioncoeffs2, approx_error_plot2

end
