"""
    rse = PartitionTreeMatrix_unbalanced_haar(matrix::Matrix{Float64},jmax::Int64;method::Symbol = :Totalvar, p::Float64 = 2.)


Perform the quadtree partition of image matrix depending on total-variation or changes on the boundary.

### Input Arguments
* `matrix`: The image matrix.
* `jmax`: The maximum level of partition.
* `method`: ':Totalvar' indicates using total variation, the other indicates using changes on the boundary.
* `p`: The power used in the criterion to balance the sizes of the sub-regions. The larger `p`, the more weight on balancing the sub-regions.

### Output Arguments
* `rse`: Region start and end. Each pair represents the start and end of a region in the matrix, corresponding to one node in the quadtree.

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function PartitionTreeMatrix_unbalanced_haar(matrix::Matrix{Float64},jmax::Int64;method::Symbol = :Totalvar, p::Float64 = 2.)
    rows,cols = size(matrix)
    N = rows*cols
    rse = Vector{Vector{Int64}}(jmax)
    rse[1] = [1,N];

    j = 1
    while j<=jmax-1
        rr = 1
        rse[j+1] = zeros(4*length(rse[j]))
        for r = 1:2:(length(rse[j]) - 1)
            rs_rows,rs_cols = ind2sub((rows,cols),rse[j][r])
            re_rows,re_cols = ind2sub((rows,cols),rse[j][r+1])

            n_rows = re_rows-rs_rows+1
            n_cols = re_cols-rs_cols+1

            n1_rows,n1_cols = PartitionTreeMatrix_unbalanced_haar_subblock(matrix[rs_rows:re_rows,rs_cols:re_cols],method = method, p = p)

            if !(n1_rows==0 || n1_rows == n_rows) && !(n1_cols==0 ||n1_cols==n_cols)
                child_rows = [rs_rows, rs_rows+n1_rows-1, rs_rows+n1_rows, re_rows, rs_rows, rs_rows+n1_rows-1, rs_rows+n1_rows, re_rows];
                child_cols = [rs_cols, rs_cols+n1_cols-1, rs_cols, rs_cols+n1_cols-1, rs_cols+n1_cols, re_cols, rs_cols+n1_cols, re_cols];
                for k = 1:8
                    rse[j+1][rr+k-1] = sub2ind((rows,cols),child_rows[k],child_cols[k])
                end
                rr = rr + 8

            elseif (n1_rows == 0 || n1_rows == n_rows) && !(n1_cols == 0 || n1_cols == n_cols)
                child_rows = [rs_rows, re_rows, rs_rows, re_rows];
                child_cols = [rs_cols, rs_cols+n1_cols-1, rs_cols+n1_cols, re_cols];
                for k = 1:4
                    rse[j+1][rr+k-1] = sub2ind((rows,cols),child_rows[k],child_cols[k])
                end
                rr = rr + 4

            elseif (n1_cols == 0 || n1_cols == n_cols) && !(n1_rows == 0 || n1_rows == n_rows)
                child_rows = [rs_rows, rs_rows+n1_rows-1, rs_rows+n1_rows, re_rows];
                child_cols = [rs_cols, re_cols, rs_cols, re_cols];
                for k = 1:4
                    rse[j+1][rr+k-1] = sub2ind((rows,cols),child_rows[k],child_cols[k])
                end
                rr = rr + 4
            else
                rse[j+1][rr:rr+1] = rse[j][r:r+1]
                rr = rr + 2
            end
        end
        j = j+1
        rse[j] = rse[j][rse[j].!=0]
    end
    return rse
end



"""
     row_cut,col_cut = PartitionTreeMatrix_unbalanced_haar_subblock(matrix::Matrix{Float64};method::Symbol = :Totalvar, p::Float64 = 2.)

Partition a matrix block into four sub-blocks.

### Input Arguments
* `matrix`: The image matrix.
* `method`: ':Totalvar' indicates using total variation, the other indicates using changes on the boundary.
* `p`: The power used in the criterion to balance the sizes of the sub-regions. The larger `p`, the more weight on balancing the sub-regions.

### Output Arguments
* `row_cut`: The cut index on the row direction.
* `col_cut`: The cut index on the col direction.

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function PartitionTreeMatrix_unbalanced_haar_subblock(matrix::Matrix{Float64};method::Symbol = :Totalvar, p::Float64 = 2.)
    m,n = size(matrix)
    Entropy = Inf

    if m>1 && n>1
        for i=1:m-1
            for j=1:n-1
                if method == :Totalvar
                    Entropy_temp = total_var(matrix[1:i,1:j])/(i*j)^p + total_var(matrix[i+1:m,1:j])/((m-i)*j)^p + total_var(matrix[1:i,j+1:n])/(i*(n-j))^p + total_var(matrix[i+1:m,j+1:n])/((m-i)*(n-j))^p;
                else
                    Entropy_temp = -( vecnorm(matrix[i,:]-matrix[i+1,:],1) + vecnorm(matrix[:,j] - matrix[:,j+1],1) )/((i*j)^(-p) + ((m-i)*j)^(-p) + (i*(n-j))^(-p) + ((m-i)*(n-j))^(-p));
                end
                if Entropy_temp < Entropy
                    Entropy = Entropy_temp
                    row_cut = i
                    col_cut = j
                end
            end
        end
    elseif m == 1 && n>1
        row_cut = 0
        for j = 1:n-1
            if method == :Totalvar
                Entropy_temp = total_var(matrix[1,1:j])/j^p + total_var(matrix[1,j+1:n])/(n-j)^p;
            else
                Entropy_temp = -norm(matrix[:,j] - matrix[:,j+1],1)*(j^(-p) + (n-j)^(-p));
            end
            if Entropy_temp < Entropy
                Entropy = Entropy_temp
                col_cut = j
            end
        end
    elseif n == 1 && m>1
        col_cut = 0
        for i = 1:m-1
            if method == :Totalvar
                Entropy_temp = total_var(matrix[1:i,1])/i^p + total_var(matrix[i+1:m,1])/(m-i)^p;
            else
                Entropy_temp = -norm(matrix[i,:]-matrix[i+1,:],1)*(i^(-p) + (m-i)^(-p));
            end
            if Entropy_temp < Entropy
                Entropy = Entropy_temp
                row_cut = i
            end
        end
    else
        row_cut = 0
        col_cut = 0
    end
    return row_cut,col_cut
end

"""
     v = total_var(matrix::Matrix{Float64})

Compute the total variation of matrix

Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function total_var(matrix::Matrix{Float64})
    m,n = size(matrix)
    row_val = abs.(matrix[1:m-1,:] - matrix[2:m,:])
    col_val = abs.(matrix[:,1:n-1] - matrix[:,2:n])
    v = norm(vcat(row_val[:],col_val[:]),1)
    return v
end
