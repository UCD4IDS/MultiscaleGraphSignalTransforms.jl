using SparseArrays
using Statistics
using LinearAlgebra
using Arpack
"""
    (pm,v) = partition_fiedler(W; method, v)

Partition the vertices of a graph according to the Fiedler vector

### Input Arguments
* `W::SparseMatrixCSC{Float64,Int}`: the edge weight matrix
* `method::Symbol`: the parition method to be used (:L or :Lrw or else ...)
* `v::Vector{Float64}`: the Fiedler vector supplied if `v` is a null vector

### Output Arguments
* `pm::Vector{Int}`: a vector of 1's and -1's
* `v::Vector{Float64}`: the Fiedler vector used for graph partitioning

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
Translated and revised by Naoki Saito, Feb. 9, 2017
"""
function partition_fiedler(W::SparseMatrixCSC{Float64,Int};
                           method::Symbol = :Lrw,
                           v::Vector{Float64} =Vector{Float64}(undef, 0))
#
# Easy case: the Fiedler vector is provided
#
    if !isempty(v)
        (pm, v) = partition_fiedler_pm(v)
        if size(W, 1) == length(v)
            # recovering the algebraic connectivity a(G)
            # Note: regardless of v coming from L or Lrw, the eigenvalue
            # computation below is fine: v'*L*v = a(G) or arw(G) because
            # for L: v'*L*v = \lambda v'*v = \lambda
            # for Lrw: vrw'*L*vrw = \lambda vrw'*D*vrw = \lambda
            # Note also, the result of the quadratic form is scalar
            # mathematically, but in julia, it's 1x1 array; so need to
            # convert it to a real scalar via `val = val[1]`.
            val = v' * (diagm(0 => vec(sum(W, dims = 1))) - W) * v; val = val[1]
            pm = partition_fiedler_troubleshooting(pm, v, W, val)
        else
            pm = partition_fiedler_troubleshooting(pm)
        end
        return pm, v
    end

#
# Now, we need to compute the Fiedler vector here. Preliminaries first.
#
    N = size(W, 1)
    eigs_flag = 0

    sigma = eps()
    cutoff = 128 # this value could be changed...

    # handle the case when there are 2 nodes
    if N == 2
        pm = [1; -1]
        v = pm / sqrt(2)
        return pm, v
    end

    # Use L (unnormalized) if specified or edge weights are very small
    if method == :L || minimum(sum(W, dims = 1)) < 10^3 * sigma
        if N > cutoff # for a relative large W
            v0 = ones(N) / sqrt(N)
            try
                # MATLAB, [v,val,eigs_flag] = eigs(diag(sum(W))-W,2,sigma,opts);
                val, vtmp = eigs(sparse(diagm(0 => vec(sum(W, dims = 1)))) - W,
                                 nev = 2, sigma = sigma, v0 = v0)
            catch emsg
                @warn("Exception in eigs(L) occurred: ", emsg)
                eigs_flag = 2
            end
            if eigs_flag == 0   # no problem in `eigs` is detected.
                val, ind = findmax(val) # val is set to be a scalar.
                v = vec(vtmp[:, ind])   # This is the Fiedler vector!
                # In MATLAB, val was returned as a diagonal matrix, so we had to
                # do all these manipulations while julia, val is simply a vector.
                # [~,ind] = max(diag(val));
                # v = v(:,ind);
                # val = val(ind,ind);
            end
        end
        if N <= cutoff || eigs_flag != 0 # if W is small or eigs had a problem,
                                         # then use full svd.
            vtmp, val = svd(diagm(0 => vec(sum(W, dims = 1))) - W) # v <-> U, val <-> S
            # diagm(vec(sum(W,1))) is Matrix{Float64} while W is SparseMatrix
            # But the results of the subtraction becomes Matrix{Float64}.
            # Also be careful here! val[end-2] == val[end-1] can happen.
            # If that is the case, vtmp[:,end-2] could also be the Fiedler vec.
            v = vec(vtmp[:, end - 1]) # v is set to be a vector.
            val = val[end - 1]       # val is set to a scalar.
            # In MATLAB, `full` was necessary to convert a sparse matrix format
            # to a usual matrix format as follows
            # [v,val,~] = svd(full(diag(sum(W))-W));
        end
    elseif method == :Lrw # Otherwise, use L_rw, which is normally preferred.
        if N > cutoff           # for a relatively large W
            v0 = ones(N) / sqrt(N)
            try
                # MATLAB: [v,val,eigs_flag] = ...
                #           eigs(diag(sum(W))-W,diag(sum(W)),2,sigma,opts);
                # This is L*v = \lambda*D*v case.
                temp = sparse(diagm(0 => vec(sum(W, dims = 1))))
                val, vtmp = eigs(temp - W, temp,
                                 nev = 2, sigma = sigma, v0 = v0)
            catch emsg
                @warn("Exception in eigs(Lrw) occurred: ", emsg)
                eigs_flag = 2
            end

            if eigs_flag == 0
                val, ind = findmax(val) # val is set to be a scalar.
                v = vtmp[:, ind]
                # MATLAB: v = (full(sum(W,2)).^(-0.5)) .* v
                v = vec((sum(W, dims = 2).^(-0.5)) .* v) # This is the Fiedler vector!
            end
        end
        if N <= cutoff || eigs_flag != 0 # if W is small or eigs had a problem,
                                         # then use full svd
            colsumW = vec(sum(W, dims = 1))
            D = sparse(diagm(0 => colsumW))
            D2 = sparse(diagm(0 => colsumW.^(-0.5)))
            vtmp, val = svd(Matrix(D2 * (D - W) * D2)) # SVD of Lsym
            # MATLAB: [v,val,~] = svd(full(...
            # bsxfun(@times,bsxfun(@times,full(sum(W,2)).^(-0.5),diag(sum(W))-W),
            #  full(sum(W,1)).^(-0.5)) ) );
            #
            # SVD is ordered in the order of nonincreasing singular values,
            # so, (end-1) entry corresponds to the algebraic connectivity.
            # Also be careful here! val[end-2] == val[end-1] can happen.
            # If that is the case, vtmp[:,end-2] could also be the Fiedler vec.
            val = val[end - 1]   # val is set to be a scalar.
            v = vtmp[:, end - 1] # v is set to be a vector.
            v = vec((sum(W, dims = 2).^(-0.5)) .* v) # This is the Fiedler vector of Lrw!
        end
    else                        # Unknown method is specified.
        error("Graph partitioning method :", method, " is not recognized!")
    end

    # partition using the Fiedler vector and troubleshoot any potential issues
    pm, v = partition_fiedler_pm(v)
    pm = partition_fiedler_troubleshooting(pm, v, W, val)
    return pm, v
end # of partition_fiedler


"""
    (pm, v) = partition_fiedler_pm(v::Vector{Float64})

    Partition an input Graph based on the Fiedler vector `v`

### Input Argument
* `v::Vector{Float64}`: the input Fiedler vector

### Output Arguments
* `pm::Vector{Int}`: the partition info
* `v::Vector{Float64}`: the output Fiedler vector if requested

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
Translated and revised by Naoki Saito, Feb. 10, 2017
"""
function partition_fiedler_pm(v::Vector{Float64})
    # define a tolerance: if abs(x) < tol, then we view x = 0.
    tol = 10^3 * eps()

    # set the first nonzero element to be positive
    row = 1
    while abs(v[row]) < tol && row < length(v)
        row = row + 1
    end
    if v[row] < 0
        v = -v
    end

    # assign each point to either region 1 or region -1, and assign any zero
    # entries to the smaller region
    if sum(v .>= tol) > sum(v .<= -tol) # more (+) than (-) entries
        pm = 2*(v .>= tol) .- 1
    else
        pm = 2*(v .<= -tol) .- 1
    end

    # make sure the first point is assigned to region 1 (not -1)
    if pm[1] < 0
        pm = -pm
    end
    return pm, v
end # of partition_fiedler_pm


"""
    pm = partition_fiedler_troubleshooting(pm,v,W,val)

Troubleshoot potential issues with the partitioning

### Input Arguments
* `pm::Vector{Int}`: an input partition info (+1 or -1) of each node
* `v::Vector{Float64}`: the Fiedler vector of the graph
* `W::SparseMatrixCSC{Float64,Int}`: an input edge weight matrix
* `val::Float64`: the algebraic connectivity (i.e., ``\\lambda_2``)

### Ouput Arguments
* `pm::Vector{Int}`: a final partition info vector

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
Translated and revised by Naoki Saito, Feb. 10, 2017
"""
function partition_fiedler_troubleshooting(pm::Vector{Int},v::Vector{Float64},
                                           W::SparseMatrixCSC{Float64,Int},
                                           val::Float64)
    # Initial setup
    N = length(pm)
    tol = 10^3 * eps()

    # If pm vector indicates no partition (single sign) or ambiguities via
    # 0 entries, then we'll trouble shoot.
    if sum(pm .< 0) == 0 || sum(pm .> 0) == 0 || sum(abs.(pm)) < N
        # Case 1: an input graph is not connected (i.e., alg. conn < tol)
        if val < tol
            pm = 2 .* (abs.(v) .> tol) .- 1
            while sum(pm .< 0) == 0 || sum(pm .> 0) == 0
                tol = 10 * tol    # relax the tolerance
                pm = 2 .* (abs.(v) .> tol) .- 1
                if tol > 1
                    pm[1:Int64(ceil(N / 2))] .= 1
                    pm[(Int64(ceil(N / 2)) + 1):N] .= -1
                end
            end
        # Case 2: it is connected, but something funny happened
        else
            pm = 2 .* (v .>= mean(v)) .- 1
            if sum(abs.(pm)) < N
                # assign the near-zero points based on the values of v
                # at their neighbor nodes
                pm0 = find(pm .== 0)
                pm[pm0] = (W[pm0, :] * v .> tol) - (W[pm0, :] * v .< -tol)
                # assign any remaining zeros to the group with fewer members
                pm[find(pm .== 0)] = (sum(pm .> 0) - sum(pm .< 0)) .>= 0
            end
        end
        # if one region has no points after all the above processing
        if sum(pm .< 0) == 0 || sum(pm .> 0) == 0
            pm[1:Int64(ceil(N / 2))] .= 1
            pm[Int64(ceil(N / 2) + 1):N] .= -1
        end
    end

# make sure that the first point is assigned as a 1
    if pm[1] < 0
        pm = -pm
    end
    return pm
end # of partition_fiedler_troubleshooting

"""
    pm = partition_fiedler_troubleshooting(pm)

Troubleshoot potential issues with the partitioning

### Input Arguments
* `pm::Vector{Int}`: an input partition info (+1 or -1) of each node

### Ouput Arguments
* `pm::Vector{Int}`: a final partition info vector

Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
Translated and revised by Naoki Saito, Feb. 10, 2017
"""
function partition_fiedler_troubleshooting(pm::Vector{Int})
    N = length(pm)
    if sum(pm .< 0) == 0 || sum(pm .> 0) == 0
        pm[1:Int(ceil(N / 2))] = 1
        pm[Int(ceil(N / 2) + 1):N] = -1
    end
    # make sure that the first point is assigned as a 1
    if pm[1] < 0
        pm = -pm
    end
    return pm
end
