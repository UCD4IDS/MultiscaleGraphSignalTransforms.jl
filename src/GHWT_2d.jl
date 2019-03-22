module GHWT_2d

include("utils.jl")

using ..GraphSignal, ..GraphPartition, ..BasisSpecification, ..GHWT, LinearAlgebra, SparseArrays

include("common.jl")

export ghwt_analysis_2d, ghwt_bestbasis_2d, ghwt_synthesis_2d




"""
    dmatrix = ghwt_analysis_2d(matrix::Array{Float64,2}, GProws::GraphPart, GPcols::GraphPart)
For a matrix, equipped with row and column weight recursive partitionings
and weight matrices, generate the redundant matrix of GHWT
expansion coefficients.

### Input Arguments
   `matrix`              the matrix to be analyzed
   `GProws`              the recursive partitioning on the rows
   `GPcols`              the recursive partitioning on the columns

### Output Arguments
   `dmatrix`                the GHWT dictionary expansion coefficients

Copyright 2019 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function ghwt_analysis_2d(matrix::Array{Float64,2}, GProws::GraphPart, GPcols::GraphPart)
    if isempty(GProws.tag)
        ghwt_core!(GProws)
    end

    if isempty(GPcols.tag)
        ghwt_core!(GPcols)
    end

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

    return Array{Float64,2}(dmatrix)
end


"""
    dvec, BSrows, BScols = ghwt_2d_bestbasis(matrix::Array{Float64,2},GProws::GraphPart,GPcols::GraphPart,
    costfun_rows::Any = 1.0, costfun_cols::Any = 1.0, flatten_rows::Any = 1.0, flatten_cols::Any = 1.0)

For a matrix, equipped with row and column weight recursive partitionings
and weight matrices, generate the (non-redundant) matrix of GHWT
expansion coefficients, using the best basis algorithm to select best bases for the rows and columns

### Input Argument
*   `matrix`              the matrix to be analyzed
*   `GProws`              the recursive partitioning on the rows
*   `GPcols`              the recursive partitioning on the columns
*   `costfun_rows`        the cost functional to be used for the rows
*   `costfun_cols`        the cost functional to be used for the columns
*   `flatten_rows`        the method for flattening vector-valued data to scalar-valued data for the rows
*   `flatten_cols`        the method for flattening vector-valued data to scalar-valued data for the columns

### Output Argument
*   `dvec`                the GHWT expansion coefficients (not redundant)
*   `BSrows`              the best basis on the rows
*   `BScols`              the best basis on the columns


Copyright 2019 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""
function ghwt_bestbasis_2d(matrix::Array{Float64,2},GProws::GraphPart,GPcols::GraphPart;
    costfun_rows::Any = 1.0, costfun_cols::Any = 1.0, flatten_rows::Any = 1.0, flatten_cols::Any = 1.0)


    # generate GraphSig objects for the matrix
    rows,cols = size(matrix)
    Grows = GraphSig(spzeros(rows,rows),f = matrix)
    Gcols = GraphSig(spzeros(cols,cols),f = Matrix{Float64}(matrix'))

    # analyze the data matrix using the rows
    dmatrix = ghwt_analysis!(Grows, GP = GProws)

    # analyze the data matrix using the columns
    dmatrix2 = ghwt_analysis!(Gcols, GP = GPcols)

    # find the row and column best bases
    dvec,BSrows = ghwt_bestbasis(dmatrix,GProws,cfspec = costfun_rows,flatten = flatten_rows)
    ~,BScols = ghwt_bestbasis(dmatrix2,GPcols,cfspec = costfun_cols,flatten = flatten_cols)

    # analyze the row best-basis coefficient matrix using the columns
    Gcols.f = Matrix{Float64}(dvec')
    dmatrix = ghwt_analysis!(Gcols, GP = GPcols)

    # extract the col best-basis coefficients
    dvec = Matrix{Float64}(dmatrix2dvec(dmatrix, GPcols, BScols)')

    return dvec, BSrows, BScols

end



"""
matrix = ghwt_synthesis_2d(dvec::Matrix{Float64},GProws::GraphPart,GPcols::GraphPart,BSrows::BasisSpec,BScols::BasisSpec)

Given a vector of GHWT expansion coefficients and info about the graph
partitioning and the choice of basis, reconstruct the signal

### Input Arguments
*   `dvec`        the expansion coefficients corresponding to the chosen basis
*   `GProws`      the recursive partitioning on the rows
*   `GPcols`      the recursive partitioning on the columns
*   `BSrows`      the best basis on the rows
*   `BScols`      the best basis on the columns

### Output Arguments
*   matrix      the reconstructed matrix

Copyright 2019 The Regents of the University of California

Implemented by Yiqun Shao (Adviser: Dr. Naoki Saito)
"""

function ghwt_synthesis_2d(dvec::Matrix{Float64},GProws::GraphPart,GPcols::GraphPart,BSrows::BasisSpec,BScols::BasisSpec)
    # reconstruct
    dvec = ghwt_synthesis(Matrix{Float64}(dvec'), GPcols, BScols)
    matrix = ghwt_synthesis(Matrix{Float64}(dvec'), GProws, BSrows)
    return matrix
end


end
