"""
    gplot(A, xyz; plotp = true, style = :auto, width = 2, color = :blue,
          shape = :none, mwidth = 1, mcolor = :blue, malpha = 0, grid = false)

GPLOT Plot graph, as in graph theory. GPLOT(A,xyz,...) plots the graph
specified by A and xyz.

### Input Arguments
* `A::SparseMatrixCSC{Float64,Int}`: the adjacency matrix of a graph `G`
* `xyz::Matrix{Float64}`: The coordinates array, `xyz`, is an n-by-2 or n-by-3 matrix with the position for node i in the i-th row, xyz[i,:] = [x[i] y[i]] or xyz[i,:] = [x[i] y[i] z[i]].
* `plotp::Bool`: if the plot is made (default) or return the X, Y, Z arrays
* `style::Symbol`: choose from Symbol[:auto,:solid,:dash,:dot,:dashdot]
* `width::Int`: in pixels
* `color::Symbol`: of line
* `shape::Symbol`: choose from Symbol[:none,:auto,:circle,:rect,:star5,:diamond,
                                      :hexagon,:cross,:xcross,:utriangle,
                                      :dtriangle,:pentagon,:heptagon,:octagon,
                                      :star4,:star6,:star7,:star8,:vline,:hline]
* `mwidth::Int`: radius in pixels
* `mcolor::Symbol`: of marker
* `malpha::Float64`: opacity of marker interior; choose from [0,1]
* `grid::Bool`: a flag to show grid lines (default: false)

### Output Arguments
* `X::Matrix{Float64}`: Nan-punctuated X coordinate vector
* `Y::Matrix{Float64}`: Nan-punctuated Y coordinate vector
* `Z::Matrix{Float64}`: Nan-punctuated Z coordinate vector

 (X,Y)= GPLOT(A,xyz,plotp=false,...) or (X,Y,Z) = GPLOT(A,xyz,plotp=false,...)
return the NaN-punctuated vectors X and Y or X, Y and Z without actually
generating a plot. These vectors can be used to generate the plot at a later
time with PLOT or PLOT3D if desired.

A backward-compatible elaboration of Mathworks's gplot
that uses 3D data (if available) when the plot is rotated.
Robert Piche, Tampere Univ. of Tech., 2005

Translated into julia by Naoki Saito, Dec. 21, 2016.
Note that we should still consider the keywords organized as
`function gplot(A,xyz,kw...)`

Nicholas Hausch edits:
Must install Plots package
For 3D rotation, use `plotlyjs()` backend

Revised by Naoki Saito, Feb. 17, 2017
"""
function gplot(A::SparseMatrixCSC{Float64,Int}, xyz::Matrix{Float64};
               plotp::Bool = true, style::Symbol = :auto, width::Int = 2,
               color::Symbol = :blue, shape::Symbol = :none, mwidth::Int = 1,
               mcolor::Symbol = :blue, malpha::Float64 = 0.0, grid::Bool = false)

    # Extract the node indices
    (i,j)=ind2sub(size(A), find(A .!= 0.0));
    p = sortperm(max(i,j))
    i = i[p]
    j = j[p]

    # Find out whether the graph is 2-D or 3-D
    nspacedim = size(xyz, 2)
    if nspacedim == 1
        xyz = hcat(xyz, zeros(size(xyz, 1)))
        nspacedim = 2
    end

    # Create a long, NaN-separated list of line segments
    X = [ xyz[i,1] xyz[j,1] Base.fill(NaN,size(i)) ]'
    Y = [ xyz[i,2] xyz[j,2] Base.fill(NaN,size(i)) ]'
    X = X[:]
    Y = Y[:]
    if nspacedim == 3
        Z = [ xyz[i,3] xyz[j,3] Base.fill(NaN,size(i)) ]'
        Z = Z[:]
    end

    # Finally, do the work!
    if plotp                    # plot the line-segments
        if nspacedim < 3        # 2D
            plot(X, Y, line=(style,width,color), marker=(shape,mwidth,mcolor),
                 markeralpha=malpha, markerstrokealpha=1, grid=grid)
        else                    # 3D
            plot(X, Y, Z, line=(style,width,color), marker=(shape,mwidth,mcolor),
                 markeralpha=malpha, markerstrokealpha=1, grid=grid)
        end
    else                        # Return the X, Y, Z array if desired.
        if nspacedim < 3
            return X, Y
        else
            return X, Y, Z
        end
    end
end
