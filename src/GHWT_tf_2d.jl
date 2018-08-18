include("PartitionTreeMatrixDhillon.jl")
include("partition_fiedler.jl")

"""
    tf2d_init(tag::Matrix{Int64},tag_r::Matrix{Int64})

Storing the relation between the coefficient location in the 2D coefficients matrix of all levels with the (level, tag, region) information

### Input Arguments
* `tag`: The matrix of the tag information
* `tag_r`: The matrix of the region information

### Output Arguments
* `tag2ind`: tag2ind[(j,k,l)] is the location of the column of coefficients on level j, region k and tag l
* `ind2tag`: key and item of dicionary `tag2ind` interchanged

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function tf2d_init(tag::Matrix{Int64},tag_r::Matrix{Int64})
    (m,n) = size(tag)
    tag2ind = Dict{Tuple{Int64,Int64,Int64},Int64}()
    ind2tag = Dict{Int64,Tuple{Int64,Int64,Int64}}()
    for j = 1:n
        for i = 1:m
            tag2ind[(j,tag_r[i,j],tag[i,j])] = i + (j-1)*m
            ind2tag[i+(j-1)*m] = (j,tag_r[i,j],tag[i,j])
        end
    end
    return tag2ind, ind2tag
end




"""
    dmatrix_new,tag2ind_new,ind2tag_new,tag_tf = tf_core_2d_col(dmatrix::Matrix{Float64},tag2ind::Dict{Tuple{Int64,Int64,Int64},Int64},ind2tag::Dict{Int64,Tuple{Int64,Int64,Int64}},jmax::Int64)

One forward iteration of 2-d time-frequency adapted ghwt on the column direction


### Input Arguments
* `dmatrix`: The cost-functional values of coeffcients
* `tag2ind`: tag2ind[(j,k,l)] is the location of the column of coefficients on level j, region k and tag l
* `ind2tag`: key and item of dicionary `tag2ind` interchanged


### Output Arguments
* `dmatrix_new`: The cost-functional values of coeffcients
* `tag2ind_new`: tag2ind_new[(j,k,l)] is the location of the column of coefficients on level j, region k and tag l
* `ind2tag_new`: key and item of dicionary `tag2ind` interchanged
* `tag_tf`: recording the time-or-frequency information


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function tf_core_2d_col(dmatrix::Matrix{Float64},tag2ind::Dict{Tuple{Int64,Int64,Int64},Int64},ind2tag::Dict{Int64,Tuple{Int64,Int64,Int64}},jmax::Int64)
    # power of the cost-functional
    costp = 1
    p,q = size(dmatrix)
    dmatrix_new = zeros(p,q)

    # Initialization
    tag_tf = zeros(Bool,(p,q))
    tag2ind_new = Dict{Tuple{Int64,Int64,Int64},Int64}()
    ind2tag_new = Dict{Int64,Tuple{Int64,Int64,Int64}}()

    s = 1
    # Iterate through columns
    for i = 1:q
        # Find out the corresponding level,region,tag information
        j,k,l = ind2tag[i]
        # only look at the entry with even tag l to avoid duplication
        if l%2 == 0 && j < jmax
            # search for the (j,k,l+1) entry.
            # the (j,k,l) and (j,k,l+1) span the same subspace as (j+1,2*k,l/2) and (j+1,2*k+1,l/2)
            # (j,k,l) and (j,k,l+1) are `frequency-direction`
            # (j,k,l) and (j,k,l+1) are `time-direction`
            if haskey(tag2ind,(j,k,l+1)) # check for degenerate case ((j,k,l+1) doesn't exist)
                freqcos = abs.(dmatrix[:,i]).^costp + abs.(dmatrix[:,tag2ind[(j,k,l+1)]]).^costp
            else
                freqcos = abs.(dmatrix[:,i]).^costp
            end
            if ~haskey(tag2ind,(j+1,2*k,l/2)) # check for degenerate case ((j+1,2*k,l/2) or (j+1,2*k+1,l/2) doesn't exist)
                timecos = abs.(dmatrix[:,tag2ind[(j+1,2*k+1,l/2)]]).^costp
            elseif ~haskey(tag2ind,(j+1,2*k+1,l/2))
                timecos = abs.(dmatrix[:,tag2ind[(j+1,2*k,l/2)]]).^costp
            else
                timecos = abs.(dmatrix[:,tag2ind[(j+1,2*k+1,l/2)]]).^costp + abs.(dmatrix[:,tag2ind[(j+1,2*k,l/2)]]).^costp
            end
            tag_tf[:,s] = timecos.>=freqcos
            dmatrix_new[:,s] = tag_tf[:,s].*freqcos + .~tag_tf[:,s].*timecos
            ind2tag_new[s] = (j,k,l/2)
            tag2ind_new[(j,k,l/2)] = s
            s = s + 1
        end
    end
    dmatrix_new = dmatrix_new[:,1:s-1]
    tag_tf = tag_tf[:,1:s-1]
    return dmatrix_new,tag2ind_new,ind2tag_new,Matrix{UInt8}(tag_tf)
end


"""
    dmatrix_new,tag2ind_new,ind2tag_new,tag_tf = tf_core_2d_row(dmatrix::Matrix{Float64},tag2ind::Dict{Tuple{Int64,Int64,Int64},Int64},ind2tag::Dict{Int64,Tuple{Int64,Int64,Int64}},jmax::Int64)

Almost the same as tf_core_2d_col. But all the operations on matrix are row-wise.
Due to the column-based matrix feature of Julia language. We are currently only using the tf_core_2d_col here.


### Input Arguments
* `dmatrix`: The cost-functional values of coeffcients
* `tag2ind`: tag2ind[(j,k,l)] is the location of the column of coefficients on level j, region k and tag l
* `ind2tag`: key and item of dicionary `tag2ind` interchanged


### Output Arguments
* `dmatrix_new`: The cost-functional values of coeffcients
* `tag2ind_new`: tag2ind_new[(j,k,l)] is the location of the row of coefficients on level j, region k and tag l
* `ind2tag_new`: key and item of dicionary `tag2ind` interchanged
* `tag_tf`: recording the time-or-frequency information


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function tf_core_2d_row(dmatrix::Matrix{Float64},tag2ind::Dict{Tuple{Int64,Int64,Int64},Int64},ind2tag::Dict{Int64,Tuple{Int64,Int64,Int64}},jmax::Int64)
    costp = 1
    p,q = size(dmatrix)
    dmatrix_new = zeros(p,q)
    tag_tf = zeros(Bool,(p,q))
    tag2ind_new = Dict{Tuple{Int64,Int64,Int64},Int64}()
    ind2tag_new = Dict{Int64,Tuple{Int64,Int64,Int64}}()

    s = 1
    for i = 1:p
        j,k,l = ind2tag[i]
        if l%2 == 0 && j < jmax
            if haskey(tag2ind,(j,k,l+1))
                freqcos = abs.(dmatrix[i,:]).^costp + abs.(dmatrix[tag2ind[(j,k,l+1)],:]).^costp
            else
                freqcos = abs.(dmatrix[i,:]).^costp
            end
            if ~haskey(tag2ind,(j+1,2*k,l/2))
                timecos = abs.(dmatrix[tag2ind[(j+1,2*k+1,l/2)],:]).^costp
            elseif ~haskey(tag2ind,(j+1,2*k+1,l/2))
                timecos = abs.(dmatrix[tag2ind[(j+1,2*k,l/2)],:]).^costp
            else
                timecos = abs.(dmatrix[tag2ind[(j+1,2*k+1,l/2)],:]).^costp + abs.(dmatrix[tag2ind[(j+1,2*k,l/2)],:]).^costp
            end
            tag_tf[s,:] = timecos.>=freqcos
            dmatrix_new[s,:] = tag_tf[s,:].*freqcos + .~tag_tf[s,:].*timecos
            ind2tag_new[s] = (j,k,l/2)
            tag2ind_new[(j,k,l/2)] = s
            s = s + 1
        end
    end
    dmatrix_new = dmatrix_new[1:s-1,:]
    tag_tf = tag_tf[1:s-1,:]
    return dmatrix_new,tag2ind_new,ind2tag_new,Matrix{UInt8}(tag_tf)
end



"""
    Bbasis, infovec = ghwt_tf_bestbasis_2d_new(dmatrix::Matrix{Float64},GProws::GraphPart,GPcols::GraphPart)

Implementation of 2d time-frequency adapted GHWT method.
Modified from the idea in paper "Image compression with adaptive Haar-Walsh tilings" by Maj Lindberg and Lars F. Villemoes.

### Input Arguments
* `dmatrix`: The ghwt expanding coefficients of matrix on all levels contatenated in a 2d matrix.
* `GProws`: Corresponding to the affinity matrix on rows.
* `GPcols`: Corresponding to the affinity matrix on cols.


### Output Arguments
* `Bbasis`: The coefficients of the best-basis.
* `infovec`: [infovec[i,1],infovec[i,2]] is the location of Bbasis[i] in dmatrix.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function ghwt_tf_bestbasis_2d(dmatrix::Matrix{Float64},GProws::GraphPart,GPcols::GraphPart)
    # Initialization
    tag_r_row = rs_to_region(GProws.rs, GProws.tag)
    tag_r_col = rs_to_region(GPcols.rs, GPcols.tag)
    tag2ind_row,ind2tag_row = tf2d_init(Matrix{Int64}(GProws.tag),Matrix{Int64}(tag_r_row))
    tag2ind_col,ind2tag_col = tf2d_init(Matrix{Int64}(GPcols.tag),Matrix{Int64}(tag_r_col))

    jmax_col = size(GPcols.tag,2)
    jmax_row = size(GProws.tag,2)

    # (level,region,tag) information on the columns
    TAG2IND_col = Vector(jmax_col)
    TAG2IND_col[1] = tag2ind_col
    IND2TAG_col = Vector(jmax_col)
    IND2TAG_col[1] = ind2tag_col

    # (level,region,tag) information on the rows
    TAG2IND_row = Vector(jmax_row)
    TAG2IND_row[1] = tag2ind_row
    IND2TAG_row = Vector(jmax_row)
    IND2TAG_row[1] = ind2tag_row

    # recording if time or frequency direction is chosen
    TAG_tf = Matrix(jmax_row,jmax_col)

    # storing the cost functional values
    DMATRIX = Matrix(jmax_row, jmax_col)

    DMATRIX[1,1] = dmatrix

    # comparison in the row direction
    for i=1:(jmax_row - 1)
        DMATRIX[i+1,1], TAG2IND_row[i+1], IND2TAG_row[i+1], TAG_tf[i+1,1] = tf_core_2d_col(DMATRIX[i,1]', TAG2IND_row[i], IND2TAG_row[i],jmax_row+1-i)
        TAG_tf[i+1,1] = (TAG_tf[i+1,1] + UInt8(2))'
        DMATRIX[i+1,1] = DMATRIX[i+1,1]'
    end

    # comparison in the column direction
    for i=1:(jmax_col - 1)
        DMATRIX[1,i+1], TAG2IND_col[i+1], IND2TAG_col[i+1], TAG_tf[1,i+1] = tf_core_2d_col(DMATRIX[1,i],TAG2IND_col[i], IND2TAG_col[i], jmax_col+1-i)
    end

    #DMATRIX[1,1] = 0 # release the memory
    # comparison in both directions
    for i = 2:jmax_row
        for j = 2:jmax_col
            (temp_col,_,_,temp_tf_col) = tf_core_2d_col(DMATRIX[i,j-1], TAG2IND_col[j-1], IND2TAG_col[j-1],jmax_col + 2 - j)
            (temp_row,_,_,temp_tf_row) = tf_core_2d_col(DMATRIX[i-1,j]', TAG2IND_row[i-1], IND2TAG_row[i-1],jmax_row + 2 - i)
            temp_row = temp_row'
            temp_tf_row = (temp_tf_row + UInt8(2))'

            row_or_col = temp_col.<=temp_row
            DMATRIX[i,j] = row_or_col.*temp_col + .~row_or_col.*temp_row
            TAG_tf[i,j] = Matrix{UInt8}(row_or_col).*temp_tf_col + Matrix{UInt8}(.~row_or_col).*temp_tf_row

            #DMATRIX[i-1,j] = 0
        end
        #DMATRIX[i,1] = 0 # release the memory
    end

    #for j = 2: jmax_col
        #DMATRIX[jmax_row, j] = 0 # release the memory
    #end
    #recover the bestbasis
    infovec = Array{Int}([jmax_row; jmax_col; 1; 1])
    for iter = 1:(jmax_row + jmax_col - 2)
        newinfovec = -1*ones(Int64,(4,2*size(infovec,2)))
        for h = 1:size(infovec,2)
            m = infovec[1,h]
            n = infovec[2,h]
            p = infovec[3,h]
            q = infovec[4,h]

            #if the bestbasis comes from column direction
            if TAG_tf[m,n][p,q] == 0 #time domain
                j,k,l = IND2TAG_col[n][q]
                if haskey(TAG2IND_col[n-1],(j+1,2*k,l))
                    q1 = TAG2IND_col[n-1][(j+1,2*k,l)]
                    newinfovec[:,2*h-1] = [m;n-1;p;q1]
                end
                if haskey(TAG2IND_col[n-1],(j+1,2*k+1,l))
                    q2 = TAG2IND_col[n-1][(j+1,2*k+1,l)]
                    newinfovec[:,2*h] = [m;n-1;p;q2]
                end
            end

            if TAG_tf[m,n][p,q] == 1 # frequency domain
                j,k,l = IND2TAG_col[n][q]
                if haskey(TAG2IND_col[n-1],(j,k,2*l))
                    q1 = TAG2IND_col[n-1][(j,k,2*l)]
                    newinfovec[:,2*h-1] = [m;n-1;p;q1]
                end
                if haskey(TAG2IND_col[n-1],(j,k,2*l+1))
                    q2 = TAG2IND_col[n-1][(j,k,2*l+1)]
                    newinfovec[:,2*h] = [m;n-1;p;q2]
                end
            end

            #if the bestbasis comes from row direction
            if TAG_tf[m,n][p,q] == 2 # time domain
                j,k,l = IND2TAG_row[m][p]
                if haskey(TAG2IND_row[m-1],(j+1,2*k,l))
                    p1 = TAG2IND_row[m-1][(j+1,2*k,l)]
                    newinfovec[:,2*h-1] = [m-1;n;p1;q]
                end
                if haskey(TAG2IND_row[m-1],(j+1,2*k+1,l))
                    p2 = TAG2IND_row[m-1][(j+1,2*k+1,l)]
                    newinfovec[:,2*h] = [m-1;n;p2;q]
                end
            end

            if TAG_tf[m,n][p,q] == 3 # frequency domain
                j,k,l = IND2TAG_row[m][p]
                if haskey(TAG2IND_row[m-1],(j,k,2*l))
                    p1 = TAG2IND_row[m-1][(j,k,2*l)]
                    newinfovec[:,2*h-1] = [m-1;n;p1;q]
                end
                if haskey(TAG2IND_row[m-1],(j,k,2*l+1))
                    p2 = TAG2IND_row[m-1][(j,k,2*l+1)]
                    newinfovec[:,2*h] = [m-1;n;p2;q]
                end
            end
        end
        infovec = newinfovec[:,newinfovec[1,:].!=-1]
    end
return dmatrix[sub2ind(size(dmatrix),infovec[3,:],infovec[4,:])]',infovec[3:4,:]
end




"""
    ghwt_tf_synthesis_2d_core!(dmatrix::Array{Float64,3},tag::Array{Int64,2},rs::Array{Int64,2})

Core function of ghwt_tf_synthesis_2d. Synthesize on column direction

### Input Arguments
* `dmatrix`: Same size as dmatrix, but only selected basis vectors are nonzero with expanding coeffcients.
* `tag`: Tag information in GP.
* `rs`: Rs information in GP.


### Output Arguments
* `dmatrix`: Synthesized on column direction.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function ghwt_tf_synthesis_2d_core!(dmatrix::Array{Float64,3},tag::Array{Int64,2},rs::Array{Int64,2})
    jmax = size(tag,2)
    for j = 1:jmax-1
        regioncount = countnz(rs[:,j]) - 1
        for r = 1:regioncount
            # the index that marks the start of the first subregion
            rs1 = rs[r,j]

            # the index that is one after the end of the second subregion
            rs3 = rs[r+1,j]

            #the number of points in the current region
            n = rs3 - rs1

            # proceed if not all coefficients are zero in this region
            if countnz(dmatrix[rs1:rs3-1,j,:]) > 0

                if n == 1
                    ### SCALING COEFFICIENT (n == 1)
                    dmatrix[rs1,j+1,:] = dmatrix[rs1,j,:] + dmatrix[rs1,j+1,:]
                elseif n > 1
                    # the index that marks the start of the second subregion
                    rs2 = rs1+1;
                    while rs2 < rs3 && tag[rs2,j+1] != 0
                        rs2 = rs2+1;
                    end

                    # the parent region is a copy of the subregion
                    if rs2 == rs3
                        dmatrix[rs1:rs3-1,j+1,:] = dmatrix[rs1:rs3-1,j,:] + dmatrix[rs1:rs3-1,j+1,:]

                        # the parent region has 2 child regions
                    else
                        # the number of points in the first subregion
                        n1 = rs2 - rs1

                        # the number of points in the second subregion
                        n2 = rs3 - rs2

                        ### SCALING COEFFICIENTS
                        dmatrix[rs1,j+1,:] = ( sqrt(n1)*dmatrix[rs1,j,:] + sqrt(n2)*dmatrix[rs1+1,j,:] )/sqrt(n) + dmatrix[rs1,j+1,:]
                        dmatrix[rs2,j+1,:] = ( sqrt(n2)*dmatrix[rs1,j,:] - sqrt(n1)*dmatrix[rs1+1,j,:] )/sqrt(n) + dmatrix[rs2,j+1,:]

                        # HAAR-LIKE & WALSH-LIKE COEFFICIENTS

                        # search through the remaining coefficients in each subregion
                        parent = rs1+2
                        child1 = rs1+1
                        child2 = rs2+1
                        while child1 < rs2 || child2 < rs3
                            # subregion 1 has the smaller tag
                            if child2 == rs3 || (tag[child1,j+1] < tag[child2,j+1] && child1 < rs2)
                                dmatrix[child1,j+1,:] = dmatrix[parent,j,:] + dmatrix[child1,j+1,:]
                                child1 = child1+1;
                                parent = parent+1;

                                # subregion 2 has the smaller tag
                            elseif child1 == rs2 || (tag[child2,j+1] < tag[child1,j+1] && child2 < rs3)
                                dmatrix[child2,j+1,:] = dmatrix[parent,j,:] + dmatrix[child2,j+1,:]
                                child2 = child2+1;
                                parent = parent+1;

                                # both subregions have the same tag
                            else
                                dmatrix[child1,j+1,:] = ( dmatrix[parent,j,:] + dmatrix[parent+1,j,:] )/sqrt(2) + dmatrix[child1,j+1,:]
                                dmatrix[child2,j+1,:] = ( dmatrix[parent,j,:] - dmatrix[parent+1,j,:] )/sqrt(2) + dmatrix[child2,j+1,:]
                                child1 = child1+1;
                                child2 = child2+1;
                                parent = parent+2;
                            end
                        end
                    end
                end
            end
        end
    end
end


"""
    matrix_r = ghwt_synthesis_2d(dmatrix::Matrix{Float64}, GProws::GraphPart, GPcols::GraphPart)

Synthesis the matrix from the coefficients of selected basis vectors.

### Input Arguments
* `dmatrix`: Only selected basis vectors are nonzero with expanding coeffcients.
* `GProws`: Corresponding to the affinity matrix on rows.
* `GPcols`: Corresponding to the affinity matrix on cols.


### Output Arguments
* `matrix_r`: Synthesized matrix


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function ghwt_synthesis_2d(dmatrix::Matrix{Float64}, GProws::GraphPart, GPcols::GraphPart)
    tag_col, rs_col = GPcols.tag, GPcols.rs
    tag_row, rs_row = GProws.tag, GProws.rs
    fcols, jmax_col = size(tag_col);
    frows, jmax_row = size(tag_row);
    matrix_r = dmatrix';
    matrix_r = reshape(matrix_r, fcols, jmax_col, frows*jmax_row);
    ghwt_tf_synthesis_2d_core!(matrix_r, Matrix{Int64}(tag_col), Matrix{Int64}(rs_col));
    matrix_r = matrix_r[:,end,:];
    #matrix_r = reshape(matrix_r, fcols, frows*jmax_row);
    matrix_r = reshape(matrix_r',frows,jmax_row, fcols);
    ghwt_tf_synthesis_2d_core!(matrix_r, Matrix{Int64}(tag_row), Matrix{Int64}(rs_row));
    matrix_r = matrix_r[:,end,:];
    #matrix_r = reshape(matrix_r,frows, fcols);
    return matrix_r
end



"""
    dmatrix, GProws, GPcols = ghwt_tf_init_2d(matrix::Matrix{Float64})

    Partition matrix first to get GProws and GPcols. Then expand matrix
    in two directions to get dmatrix

### Input Arguments
* `matrix::Matrix{Float64}`: an input matrix

### Output Argument
* `GProws::GraphPart`: partitioning using rows as samples
* `GPcols::GraphPart`: partitioning using cols as samples
* `dmatrix::matrix{Float64}`: expansion coefficients of matrix
"""
function ghwt_tf_init_2d(matrix::Matrix{Float64})
  #Partition matrix recuisively using Dhillon's method
  GProws, GPcols = PartitionTreeMatrixDhillon(matrix)
  ghwt_core!(GProws)
  ghwt_core!(GPcols)


  matrix_re = matrix[GProws.ind, GPcols.ind]

  (N, jmax_row) = size(GProws.rs)
  N = N-1

  # expand on the row direction
  (frows, fcols) = size(matrix_re)
  dmatrix = zeros((N,jmax_row,fcols))
  dmatrix[:,jmax_row,:] = matrix_re

  ghwt_core!(GProws, dmatrix)

  dmatrix = reshape(dmatrix, (frows*jmax_row, fcols))'

  # expand on the column direction
  (N, jmax_col) = size(GPcols.rs)
  N = N - 1
  dmatrix2 = zeros(N,jmax_col,size(dmatrix,2))
  dmatrix2[:, jmax_col, :] = dmatrix
  ghwt_core!(GPcols, dmatrix2)

  dmatrix = reshape(dmatrix2, (fcols*jmax_col, frows*jmax_row))'

  return dmatrix, GProws, GPcols
end
