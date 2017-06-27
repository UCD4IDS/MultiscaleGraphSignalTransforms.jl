module GraphPartition

using ..GraphSignal

include("partition_fiedler.jl")
include("utils.jl")

export GraphPart, partition_tree_fiedler

# GraphPart data structure and constructors
"""
    GP = GraphPart{Tl, Ts}(ind::Vector{Tl}, rs::Matrix{Tl}, tag::Matrix{Ts}, compinfo::Matrix{Tl}, rsf2c::Matrix{Tl}, tagf2c::Matrix{Ts}, compinfof2c::Matrix{Tl}, method::Vector{Symbol})

is a data structure for a GraphPart object containing the following fields:
* `ind::Vector{Tl}`: ordering of the indices on the finest level
* `rs::Matrix{Tl}`: regionstarts (coarse-to-fine) <==> the index in `ind` of the first point in region number `i` is `rs[i]`
* `tag::Matrix{Ts}`: tag info for the GHWT coarse-to-fine basis
* `compinfo::Matrix{Tl}`: indicates whether the coefficient was formed from 2 coefficents (value is nonzero) or from only 1 coefficient (value is zero); when a scaling and Haar-like coefficient are formed, their corresponding values in compinfo indicate the number of nodes in each of the 2 subregions
* `rsf2c::Matrix{Tl}`: the fine-to-coarse version of `rs`
* `tagf2c::Matrix{Ts}`: the fine-to-coarse version of `tag`
* `compinfof2c::Matrix{Tl}`: the fine-to-coarse version of `compinfo`
* `method::Symbol`: how the partition tree was constructed

The unsigned integer depends on the size of the underlying graph.

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and modified by Naoki Saito, Feb. 7, 2017
Revised for two parameters by Naoki Saito, Feb. 24, 2017
"""
type GraphPart{Tl <: Unsigned, Ts <: Unsigned}
    ind::Vector{Tl}    # ordering of the indices on the finest level
    rs::Matrix{Tl}     # `rs[i,j]` = the index in `ind` of the first
                       # point in Region `i` at level j 
    tag::Matrix{Ts}    # Ts <: Tl because max(tag) <= log2(max(T))
    compinfo::Matrix{Tl}      # coef was from 2 coefs or from 1 coef
    rsf2c::Matrix{Tl}         # f2c version of `rs`
    tagf2c::Matrix{Ts}        # f2c version of `tag`
    compinfof2c::Matrix{Tl}   # f2c version of `compinfo`
    method::Symbol           # specification of graph partition method

    # An inner constructor here.
    function GraphPart{Tl,Ts}(ind::Vector{Tl}, rs::Matrix{Tl};
                       tag::Matrix{Ts} = Matrix{Ts}(0, 0),
                       compinfo::Matrix{Tl} = Matrix{Tl}(0, 0),
                       rsf2c::Matrix{Tl} = Matrix{Tl}(0, 0),
                       tagf2c::Matrix{Ts} = Matrix{Ts}(0, 0),
                       compinfof2c::Matrix{Tl} = Matrix{Tl}(0, 0),
                       method::Symbol = :unspecified)
        
        # Sanity checks
        if Base.size(rs, 1) != Base.length(ind) + 1
            warn("size(rs,1) must be length(ind) + 1")
            warn("both rs and ind now become null arrays!")
            ind = Matrix{Tl}()
            rs = Matrix{Tl}()
        end
        if Base.length(tag) != 0 && (Base.size(rs, 1) != Base.size(tag, 1) + 1
            || Base.size(rs, 2) != Base.size(tag, 2))
            warn("tag size is inconsistent with rs size.")
            warn("tag now becomes a null array!")
            tag = Matrix{Ts}()
        end
        if Base.length(compinfo) != 0 &&
            ( Base.size(rs, 1) != Base.size(compinfo, 1) + 1
             || Base.size(rs, 2) != Base.size(compinfo, 2) )
            warn("compinfo size is inconsistent with rs size.")
            warn("compinfo now becomes a null array!")
            compinfo = Matrix{Tl}()
        end
        if Base.length(rsf2c) != 0 && Base.length(rsf2c) != Base.length(rs)
            warn("length(rsf2c) must be length(rs)")
            warn("rsf2c now becomes a null array!")
            rsf2c = Matrix{Tl}()
        end
        if Base.length(tagf2c) != 0 && Base.length(tagf2c) != Base.length(tag)
            warn("length(tagf2c) must be length(tag)")
            warn("tagf2c now becomes a null array!")
            tagf2c = Matrix{Ts}()
        end
        if Base.length(compinfof2c) != 0 &&
            Base.length(compinfof2c) != Base.length(compinfo)
            warn("length(compinfof2c) must be length(compinfo)")
            warn("compf2c now becomes a null array!")
            compf2c = Matrix{Tl}()
        end
        new(ind, rs, tag, compinfo, rsf2c, tagf2c, compinfof2c, method)
    end # of an inner constructor GraphPart

 end # of type GraphPart


"""
    GP = partition_tree_fiedler(G::GraphSig, method::Symbol = :Lrw)

 Generate a partition tree for a graph using the Fiedler vector of either
 L (the unnormalized Laplacian) or L_rw (the random-walk normalized
 Laplacian).

### Input Arguments
* `G::GraphSig`: an input GraphSig object
* `method::Symbol`: how the partition tree was constructed (default: :Lrw)

### Output Argument
* `GP::GraphPart`: an ouput GraphPart object

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated and modified by Naoki Saito, Feb. 7, 2017
"""
function partition_tree_fiedler(G::GraphSignal.GraphSig, method::Symbol = :Lrw)
    #
    # 0. Preliminary stuff
    #
    # constants
    N = G.length
    jmax = max(3 * floor(Int, log2(N)), 4) # jmax >= 4 is guaranteed.
    # This jmax is the upper bound of the true jmax; the true jmax
    # is computed as the number of columns of the matrix `rs` in the end.

    # TRACKING VARIABLES

    # `ind` records the way in which the nodes are indexed on each level
    Tl = ind_class(N)
    ind = Vector{Tl}(1:N)

    # `rs` stands for regionstarts, meaning that the index in `ind` of the first
    # point in region number `i` is `rs[i]`
    rs = zeros(Tl, N + 1, jmax)
    rs[1, :] = 1
    rs[2, 1] = N + 1

    #
    # 1. Partition the graph to yield rs, ind, and jmax
    #
    j = 1                       # j=1 in julia <=> j=0 in theory
    regioncount = 0
    rs0 = 0                     # define here for the whole loops,
                                # which differs from MATLAB.
    while regioncount < N
        regioncount = countnz(rs[:, j]) - 1 # the number of regions on level j
        if j == jmax  # add a column to rs for level j+1, if necessary
            rs = hcat(rs, vcat(1, zeros(N))) 
            jmax = jmax + 1
        end
        # for tracking the child regions
        rr = 1
        for r = 1:regioncount   # cycle through the parent regions
            rs1 = rs[r, j]      # the start of the parent region
            rs2 = rs[r + 1, j] # 1 node after the end of the parent region
            n = rs2 - rs1   # the number of nodes in the parent region
            
            if n > 1            # regions with 2 or more nodes
                indrs = ind[rs1:(rs2 - 1)]
                # partition the current region
                (pm, ) = partition_fiedler(G.W[indrs,indrs], method = method)
                # determine the number of points in child region 1
                n1 = sum(pm .> 0)
                # switch regions 1 and 2, if necessary, based on the sum of
                # edge weights to the previous child region (i.e., cousin)
                if r > 1
                    # MATLAB: sum(sum(...)) instead of sum(...)
                    if sum(G.W[ind[rs0:(rs1 - 1)], indrs[pm .> 0]]) <
                        sum(G.W[ind[rs0:(rs1 - 1)], indrs[pm .< 0]])
                        pm = -pm
                        n1 = n - n1
                    end
                end
                # update the indexing
                ind[rs1:(rs1 + n1 - 1)] = indrs[pm .> 0]
                ind[(rs1 + n1):(rs2 - 1)] = indrs[pm .< 0]
                # update the region tracking
                rs[rr + 1, j + 1] = rs1 + n1
                rs[rr + 2, j + 1] = rs2
                rr = rr + 2
                rs0 = rs1 + n1
            elseif n == 1       # regions with 1 node
                rs[rr + 1, j + 1] = rs2
                rr = rr + 1
                rs0 = rs1
            end # of if n > 1 ... elseif n==1 construct
        end # of for r=1:regioncount
        j = j + 1
    end # of while regioncount < N statement

    #
    # 2. Postprocessing
    #
    # get rid of excess columns in rs
    rs = rs[:, 1:(j - 1)]         # in MATLAB, it was rs(:,j:end) = [];
    # create a GraphPart object
    Ts = tag_class(Base.size(rs,2))
    return GraphPart{Tl, Ts}(ind, rs, method = method)
end # of function partition_tree_fiedler

end # of module GraphPartition
