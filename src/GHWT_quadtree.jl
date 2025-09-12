"""
    ind,rse,IND = PartitionTreeMatrixDhillon_QuadTree(matrix,jmax)


Perform the quadtree partition of matrix.

### Input Arguments
* `matrix`: The matrix to be partitioned
* `jmax`: Specify the max level of partitioning the matrix

### Output Arguments
* `ind`: Vector with the same length as length(matrix[:]), indicate the order of matrix after partition.
* `rse`: Region start and end. Each pair represents the start and end of a region in the matrix, corresponding to one node in the quadtree
* `IND`: IND[i] indicate the order of matrix after iteration i.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function PartitionTreeMatrixDhillon_QuadTree(matrix,jmax)

# 0.Preliminary

# Constants
rows,cols = size(matrix);
N = rows*cols
ind = Matrix{Int64}(reshape(1:N,rows,cols))
IND = zeros(rows,cols,jmax);
IND[:,:,1] = ind;

# region starts & ends
rse = Vector{Vector{Int64}}(jmax)
rse[1] = [1,N];

# 1. Partition the graph to yield rse, ind and jmax

j=1;
regioncount = 0;

while j<=jmax-1
    # add a column for j+1, if necessary

    rr = 1#for tracking the child regions

    rse[j+1] = zeros(Int64,4*length(rse[j]))

    for r = 1:2:(length(rse[j]) - 1)
        rs_rows,rs_cols = ind2sub((rows,cols),rse[j][r])
        re_rows,re_cols = ind2sub((rows,cols),rse[j][r+1])

        n_rows = re_rows-rs_rows+1
        n_cols = re_cols-rs_cols+1

        indrs = ind[rs_rows:re_rows,rs_cols:re_cols]

        if n_rows > 1 && n_cols > 1


            # compute the left and right singular vectors
            u_rows,v_cols = second_largest_singular_vectors(reshape(matrix[indrs[:]],n_rows,n_cols))

            # partition the current region
            pm,_ = partition_fiedler(SparseMatrixCSC{Float64,Int}(ones(1,1)),v = vcat(u_rows, v_cols))
            pm_rows = pm[1:n_rows]
            pm_cols = pm[n_rows+1:end]

            # determine the number of points in child region 1
            n1_rows = sum(pm_rows .> 0);
            n1_cols = sum(pm_cols .> 0);

            # update the indexing
            ind[rs_rows:re_rows,rs_cols:re_cols] = indrs[vcat(find(pm_rows .> 0),find(pm_rows.<=0)),vcat(find(pm_cols.>0),find(pm_cols.<=0))]

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

        elseif n_rows == 1 && n_cols > 1 #no bisection on rows
            pm_cols,_ = partition_fiedler(SparseMatrixCSC{Float64,Int}(ones(1,1)),v = matrix[indrs[:]])
            n1_cols = sum(pm_cols .> 0);
            ind[rs_rows:re_rows,rs_cols:re_cols] = indrs[:,vcat(find(pm_cols.>0),find(pm_cols.<=0))]

            child_rows = [rs_rows, re_rows, rs_rows, re_rows];
            child_cols = [rs_cols, rs_cols+n1_cols-1, rs_cols+n1_cols, re_cols];
            for k = 1:4
                rse[j+1][rr+k-1] = sub2ind((rows,cols),child_rows[k],child_cols[k])
            end
            rr = rr + 4


        elseif n_rows > 1 && n_cols == 1 #no bisection on cols
            pm_rows,_ = partitio_fiedler(SparseMatrixCSC{Float64,Int}(ones(1,1)),v = matrix[indrs[:]])
            n1_rows = sum(pm_rows .> 0);
            ind[rs_rows:re_rows,rs_cols:re_cols] = indrs[vcat(find(pm_rows .> 0),find(pm_rows.<=0)),:]

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
    IND[:,:,j] = ind
end
return ind,rse,IND
end





"""
    u,v = second_largest_singular_vectors(A::Matrix{Float64})


Compute the second largest left and right singular vectors of matrix.

### Input Arguments
* `A`: The matrix of which the singular vectors will be computed.

### Output Arguments
* `u`: Second largest left singular vector.
* `v`: Second largest right singular vector.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function second_largest_singular_vectors(A::Matrix{Float64})
#Compute the second largest left and right singular vectors of An = D1^-0.5*A*D2^-0.5

rows,cols = size(A)

# compute D1 and D2 and make sure there are no zeros
D1 = sum(A,2)
D1[D1 .== 0] = maximum([0.01, minimum(D1[D1 .> 0])/10])
D2 = sum(A,1)
D2[D2 .== 0] = maximum([0.01, minimum(D2[D2 .> 0])/10])

u,_,v = svd(Diagonal(D1[:].^(-0.5))*A*Diagonal(D2[:].^(-0.5)))

# extract the 2nd singular vectors and multiply by D_i^-0.5
u = Diagonal(D1[:].^-0.5)*u[:,2]
v = Diagonal(D2[:].^-0.5)*v[:,2]
return u,v
end





"""
    BBmatrix, BBmatrix_ind, tag, dmatrix = GHWT_Matrix2D_Analysis_BestBasis(matrix::Matrix{Float64},ind::Matrix{Int64} ,rse::Vector{Array{Int64,1}})


Compute the expanding coeffcients of all levels and search for the best-basis given the partition quadtree of matrix.

### Input Arguments
* `matrix`: The matrix which will be analyzed.
* `ind`: Vector with the same length as length(matrix[:]), indicate the order of matrix after partition.
* `rse`: Region start and end. Each pair represents the start and end of a region in the matrix, corresponding to one node in the quadtree

### Output Arguments
* `BBmatrix`: Coefficients of the best-basis.
* `BBmatrix_ind`: Location of the best-basis in dmatrix.
* `tag`: Tag information of the expanding coefficients.
* `dmatrix`: The expanding coefficients on all levels.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function GHWT_Matrix2D_Analysis_BestBasis(matrix::Matrix{Float64},ind::Matrix{Int64} ,rse::Vector{Array{Int64,1}})
    # cost functional
    costfun = function(x) norm(x[:],1) end

    # constants
    rows,cols = size(matrix);
    #tol = 10^4*eps()

    # 1. Recursively partition the rows and columns of the matrix in a coupled manner
    #ind,rse,_ = PartitionTreeMatrixDhillon_QuadTree_ys(matrix)
    jmax = length(rse)

    # 2. Analyze the matrix and perform the best basis search
    dmatrix = zeros(rows,cols,jmax)
    dmatrix[:,:,jmax] = matrix[ind]
    tag = zeros(Int64,rows,cols,jmax)

    for r = 1:2:length(rse[jmax])
        if rse[jmax][r] != rse[jmax][r+1]
            rs_rows,rs_cols = ind2sub((rows,cols),rse[jmax][r])
            re_rows,re_cols = ind2sub((rows,cols),rse[jmax][r+1])

            tag[rs_rows:re_rows,rs_cols:re_cols,jmax] = reshape(0:(re_rows-rs_rows+1)*(re_cols-rs_cols+1)-1,re_rows-rs_rows+1,re_cols-rs_cols+1)
            average = sum(dmatrix[rs_rows:re_rows,rs_cols:re_cols,jmax])/sqrt((re_rows-rs_rows+1)*(re_cols-rs_cols+1))
            dmatrix[rs_rows:re_rows,rs_cols:re_cols,jmax] = 0
            dmatrix[rs_rows,rs_cols,jmax] = average
        end
    end
    #best basis
    BBmatrix = dmatrix[:,:,jmax]
    BBmatrix_ind = jmax*ones(Int64,rows,cols)

    for j = jmax-1:-1:1
        rr = 1;
        for r = 1:2:length(rse[j])
            rs_rows,rs_cols = ind2sub((rows,cols),rse[j][r])
            re_rows,re_cols = ind2sub((rows,cols),rse[j][r+1])

            # the row indices being considered
            indrs_rows = rs_rows:re_rows

            # the column indices being considered
            indrs_cols = rs_cols:re_cols

            n_rows = re_rows-rs_rows+1
            n_cols = re_cols-rs_cols+1

            #copy of children
            if rse[j][r+1] == rse[j+1][rr+1]
                dmatrix[indrs_rows,indrs_cols,j] = dmatrix[indrs_rows,indrs_cols,j+1]
                tag[indrs_rows,indrs_cols,j] = tag[indrs_rows,indrs_cols,j+1]
                rr = rr+2

                #copy of two childern
            elseif rse[j][r+1] == rse[j+1][rr+3]
                a,b = zeros(Int64,4), zeros(Int64,4)
                for k = 1:4
                    a[k],b[k] = ind2sub((rows,cols),rse[j+1][rr+k-1])
                end

                region1 = dmatrix[a[1]:a[2],b[1]:b[2],j+1]
                region2 = dmatrix[a[3]:a[4],b[3]:b[4],j+1]
                tag1 = tag[a[1]:a[2],b[1]:b[2],j+1]
                tag2 = tag[a[3]:a[4],b[3]:b[4],j+1]
                region_vec, tag_vec = twoblockcombine(region1[:],region2[:],tag1[:],tag2[:]);
                tag[indrs_rows,indrs_cols,j] = reshape(tag_vec,n_rows,n_cols);
                dmatrix[indrs_rows,indrs_cols,j] = reshape(region_vec,n_rows,n_cols)
                rr = rr + 4

            elseif rse[j][r+1] == rse[j+1][rr+7]
                a,b = zeros(Int64,8), zeros(Int64,8)
                for k = 1:8
                    a[k],b[k] = ind2sub((rows,cols),rse[j+1][rr+k-1])
                end

                region11 = dmatrix[a[1]:a[2],b[1]:b[2],j+1]
                region21 = dmatrix[a[3]:a[4],b[3]:b[4],j+1]
                tag11 = tag[a[1]:a[2],b[1]:b[2],j+1]
                tag21 = tag[a[3]:a[4],b[3]:b[4],j+1]
                region12 = dmatrix[a[5]:a[6],b[5]:b[6],j+1]
                region22 = dmatrix[a[7]:a[8],b[7]:b[8],j+1]
                tag12 = tag[a[5]:a[6],b[5]:b[6],j+1]
                tag22 = tag[a[7]:a[8],b[7]:b[8],j+1]

                region_vec_1, tag_vec_1 = twoblockcombine(region11[:],region21[:],tag11[:],tag21[:])
                region_vec_2, tag_vec_2 = twoblockcombine(region12[:],region22[:],tag12[:],tag22[:])
                region_vec, tag_vec = twoblockcombine(region_vec_1,region_vec_2,tag_vec_1,tag_vec_2)

                tag[indrs_rows,indrs_cols,j] = reshape(tag_vec, n_rows,n_cols)
                dmatrix[indrs_rows,indrs_cols,j] = reshape(region_vec,n_rows,n_cols)
                rr = rr + 8

            end
            #Obtaining BBmatrix
            if costfun(BBmatrix[indrs_rows,indrs_cols]) > costfun(dmatrix[indrs_rows,indrs_cols,j])
                BBmatrix[indrs_rows,indrs_cols] = dmatrix[indrs_rows,indrs_cols,j]
                BBmatrix_ind[indrs_rows,indrs_cols] = j;

            end
        end
    end
    return BBmatrix, BBmatrix_ind, tag, dmatrix
end





"""
    region, tag = twoblockcombine(region1::Vector{Float64}, region2::Vector{Float64}, tag1::Vector{Int64}, tag2::Vector{Int64})


Compute the coefficients of coarser basis-vectors after combining finer basis-vectors supported on two sub-regions.

### Input Arguments
* `region1`: The coeffcient of finer basis-vectors supported on region 1.
* `tag1`: The tag information of finer basis-vectors supported on region 1.
* `region2`: The coeffcient of finer basis-vectors supported on region 2.
* `tag2`: The tag information of finer basis-vectors supported on region 2.

### Output Arguments
* `region`: The coeffcient of coarser basis-vectors supported on region 1 and region 2.
* `tag`: The tag information of coarser basis-vectors supported on region 1 and region 2.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function twoblockcombine(region1::Vector{Float64}, region2::Vector{Float64}, tag1::Vector{Int64}, tag2::Vector{Int64})
    # the number of points in the first subregion
    n1 = length(region1)
    # the number of points in the second subregion
    n2 = length(region2)

    n = n1+n2
    region = zeros(n)
    tag = zeros(Int64,n)

    ### SCALING COEFFICIENT (n > 1)
    region[1] = ( sqrt(n1)*region1[1] + sqrt(n2)*region2[1] )/sqrt(n)
    tag[1] = 0

    ### HAAR-LIKE COEFFICIENT
    region[2] = ( sqrt(n2)*region1[1] - sqrt(n1)*region2[1] )/sqrt(n);
    tag[2] = 1

    ### WALSH-LIKE COEFFICIENTS
    # sweep through the coefficients in subregion 1 and subregion 2

    # the index of the new coefficient(s) to be created on level j
    parent = 3

    # the index of the current coefficients in subregions 1 and 2
    child1 = 2
    child2 = 2

    while parent <= n
        # no matching coefficient (use subregion 1)
        if child1 <= n1 && (child2 == n2+1 || tag1[child1] < tag2[child2])
            region[parent] = region1[child1]
            tag[parent] = 2*tag1[child1]
            child1 = child1+1;
            parent = parent+1;

            # no matching coefficient (use subregion 2)
        elseif child2 <=n2 && (child1 == n1+1 || tag2[child2] < tag1[child1])
            region[parent] = region2[child2]
            tag[parent] = 2*tag2[child2]
            child2 = child2+1
            parent = parent+1

            # matching coefficients
        else
            region[parent] = ( region1[child1] + region2[child2] )/sqrt(2)
            tag[parent] = 2*tag[child1]

            region[parent+1] = ( region1[child1] - region2[child2] )/sqrt(2)
            tag[parent+1] = 2*tag[child1]+1

            child1 = child1+1
            child2 = child2+1
            parent = parent+2
        end
    end
    return region, tag
end





"""
    matrix_r = GHWT_Matrix2D_Analysis_Synthesis(BBmatrix::Matrix{Float64},BBmatrix_ind::Matrix{Int64},tag::Array{Int64,3} ,rse::Vector{Array{Int64,1}})


Synthesize the matrix given coefficients of selected basis-vectors.

### Input Arguments
* `BBmatrix`: Coefficients of the best-basis.
* `BBmatrix_ind`: Location of the best-basis in dmatrix.
* `tag`: Tag information of the expanding coefficients.
* `rse`: Region start and end. Each pair represents the start and end of a region in the matrix, corresponding to one node in the quadtree.

### Output Arguments
* `matrix_r`: The synthesized matrix.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function GHWT_Matrix2D_Analysis_Synthesis(BBmatrix::Matrix{Float64},BBmatrix_ind::Matrix{Int64},tag::Array{Int64,3} ,rse::Vector{Array{Int64,1}})
    costfun = function(x) norm(x[:],1) end
    # constants
    rows,cols = size(BBmatrix);
    #tol = 10^4*eps()
    jmax = length(rse)

    cmatrix = zeros(rows,cols,jmax);
    for i  = 1:rows
        for j = 1:cols
            cmatrix[i,j,BBmatrix_ind[i,j]] = BBmatrix[i,j];
        end
    end

    for j = 1:1:jmax-1
        rr = 1
        for r = 1:2:length(rse[j])
            rs_rows,rs_cols = ind2sub((rows,cols),rse[j][r])
            re_rows,re_cols = ind2sub((rows,cols),rse[j][r+1])

            if rse[j][r+1] == rse[j+1][rr+1]
                if countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j+1]) == 0 && countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]) > 0
                    cmatrix[rs_rows:re_rows,rs_cols:re_cols,j+1] = cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]
                end
                rr = rr + 2

            elseif rse[j][r+1] == rse[j+1][rr+3]
                if countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j+1]) == 0 && countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]) > 0
                    a,b = zeros(Int64,4), zeros(Int64,4)
                    for k = 1:4
                        a[k],b[k] = ind2sub((rows,cols),rse[j+1][rr+k-1])
                    end

                    region = cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]
                    tag1 = tag[a[1]:a[2],b[1]:b[2],j+1]
                    tag2 = tag[a[3]:a[4],b[3]:b[4],j+1]
                    region1_vec, region2_vec = twoblocksplit(region[:],tag1[:],tag2[:]);
                    cmatrix[a[1]:a[2],b[1]:b[2],j+1] = reshape(region1_vec,a[2]-a[1]+1,b[2]-b[1]+1)
                    cmatrix[a[3]:a[4],b[3]:b[4],j+1] = reshape(region2_vec,a[4]-a[3]+1,b[4]-b[3]+1)
                end
                rr = rr + 4
            else
                if countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j+1]) == 0 && countnz(cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]) > 0
                    a,b = zeros(Int64,8), zeros(Int64,8)
                    for k = 1:8
                        a[k],b[k] = ind2sub((rows,cols),rse[j+1][rr+k-1])
                    end
                    tag11 = tag[a[1]:a[2],b[1]:b[2],j+1]
                    tag21 = tag[a[3]:a[4],b[3]:b[4],j+1]
                    tag12 = tag[a[5]:a[6],b[5]:b[6],j+1]
                    tag22 = tag[a[7]:a[8],b[7]:b[8],j+1]
                    tag11 = tag11[:]
                    tag21 = tag21[:]
                    tag12 = tag12[:]
                    tag22 = tag22[:]
                    tag1 = twoblockcombine_tag(tag11,tag21)
                    tag2 = twoblockcombine_tag(tag12,tag22)
                    region = cmatrix[rs_rows:re_rows,rs_cols:re_cols,j]
                    region1_vec,region2_vec = twoblocksplit(region[:], tag1,tag2)
                    region11,region21 = twoblocksplit(region1_vec,tag11,tag21)
                    region12,region22 = twoblocksplit(region2_vec,tag12,tag22)
                    cmatrix[a[1]:a[2],b[1]:b[2],j+1] = reshape(region11,a[2]-a[1]+1,b[2]-b[1]+1)
                    cmatrix[a[3]:a[4],b[3]:b[4],j+1] = reshape(region21,a[4]-a[3]+1,b[4]-b[3]+1)
                    cmatrix[a[5]:a[6],b[5]:b[6],j+1] = reshape(region12,a[6]-a[5]+1,b[6]-b[5]+1)
                    cmatrix[a[7]:a[8],b[7]:b[8],j+1] = reshape(region22,a[8]-a[7]+1,b[8]-b[7]+1)
                end
                rr = rr + 8
            end
        end
    end

    matrix_r = cmatrix[:,:,end]
    for r = 1:2:length(rse[jmax])
        if rse[jmax][r] != rse[jmax][r+1]
            rs_rows,rs_cols = ind2sub((rows,cols),rse[j][r])
            re_rows,re_cols = ind2sub((rows,cols),rse[j][r+1])

            matrix_r[rs_rows:re_rows,rs_cols:re_cols] = matrix_r[rs_rows,rs_cols]/sqrt((re_rows-rs_rows+1)*(re_cols-rs_cols+1))
        end
    end
    return matrix_r
end





"""
    tag = twoblockcombine_tag(tag1::Vector{Int64},tag2::Vector{Int64})


Compute the tag information of coarser basis-vectors after combining finer basis-vectors supported on two sub-regions.

### Input Arguments
* `tag1`: The tag information of finer basis-vectors supported on region 1.
* `tag2`: The tag information of finer basis-vectors supported on region 2.

### Output Arguments
* `tag`: The tag information of coarser basis-vectors supported on region 1 and region 2.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function twoblockcombine_tag(tag1::Vector{Int64},tag2::Vector{Int64})
    # the number of points in the first subregion
    n1 = length(tag1)
    # the number of points in the second subregion
    n2 = length(tag2)

    n = n1+n2
    tag = zeros(Int64,n)

    ### SCALING COEFFICIENT (n > 1)
    tag[1] = 0

    ### HAAR-LIKE COEFFICIENT
    tag[2] = 1

    ### WALSH-LIKE COEFFICIENTS
    # sweep through the coefficients in subregion 1 and subregion 2

    # the index of the new coefficient(s) to be created on level j
    parent = 3

    # the index of the current coefficients in subregions 1 and 2
    child1 = 2
    child2 = 2

    while parent <= n
        # no matching coefficient (use subregion 1)
        if child1 <= n1 && (child2 == n2+1 || tag1[child1] < tag2[child2])
            tag[parent] = 2*tag1[child1]
            child1 = child1+1
            parent = parent+1

            # no matching coefficient (use subregion 2)
        elseif child2 <=n2 && (child1 == n1+1 || tag2[child2] < tag1[child1])
            tag[parent] = 2*tag2[child2]
            child2 = child2+1
            parent = parent+1

            # matching coefficients
        else
            tag[parent] = 2*tag[child1]
            tag[parent+1] = 2*tag[child1]+1

            child1 = child1+1
            child2 = child2+1
            parent = parent+2
        end
    end
    return tag
end




"""
    region1,region2 = twoblocksplit(region::Matrix{Float64}, tag1::Vector{Int64}, tag2::Vector{Int64})


Compute the coefficients of finer basis-vectors supported on two sub-regions given the coeffcients of coarser basis-vectors supported on both regions.

### Input Arguments
* `region`: The coeffcient of coarser basis-vectors supported on region 1 and region 2.
* `tag1`: The tag information of finer basis-vectors supported on region 1.
* `tag2`: The tag information of finer basis-vectors supported on region 2.

### Output Arguments
* `region1`: The coeffcient of finer basis-vectors supported on region 1.
* `region2`: The coeffcient of finer basis-vectors supported on region 2.


Copyright 2018 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function twoblocksplit(region::Matrix{Float64}, tag1::Vector{Int64}, tag2::Vector{Int64})
    # the number of points in the first subregion
    n1 = length(tag1)
    # the number of points in the second subregion
    n2 = length(tag2)

    n = n1+n2
    region1 = zeros(n1)
    region2 = zeros(n2);

    ### SCALING COEFFICIENT (n > 1)
    region1[1] = ( sqrt(n1)*region[1] + sqrt(n2)*region[2] )/sqrt(n);

    ### HAAR-LIKE COEFFICIENT
    region2[1] = ( sqrt(n2)*region[1] - sqrt(n1)*region[2] )/sqrt(n);


    ### WALSH-LIKE COEFFICIENTS
    # sweep through the coefficients in subregion 1 and subregion 2

    # the index of the new coefficient(s) to be created on level j
    parent = 3

    # the index of the current coefficients in subregions 1 and 2
    child1 = 2
    child2 = 2

    while parent <= n
        # no matching coefficient (use subregion 1)
        if child1 <= n1 && (child2 == n2+1 || tag1[child1] < tag2[child2])
            region1[child1] = region[parent]
            child1 = child1+1
            parent = parent+1

            # no matching coefficient (use subregion 2)
        elseif child2 <=n2 && (child1 == n1+1 || tag2[child2] < tag1[child1])
            region2[child2] = region[parent]
            child2 = child2+1
            parent = parent+1

            # matching coefficients
        else
            region1[child1] = ( region[parent] + region[parent+1] )/sqrt(2)
            region2[child2] = ( region[parent] - region[parent+1] )/sqrt(2)
            child1 = child1+1
            child2 = child2+1
            parent = parent+2
        end

    end
    return region1,region2
end
