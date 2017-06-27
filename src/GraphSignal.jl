module GraphSignal

using Distances                 # This is necessary in Ggrid
using Distributions             # This is necessary in AddNoise

export GraphSig, replace_data!, snr, gpath, ggrid, Adj2InvEuc, AddNoise

# GraphSig data structure and constructors
"""
    G = GraphSig(W, xy, f, name, plotspecs, length, dim, size)

is a data structure for a GraphSig object containing the following fields:
* `W::SparseMatrixCSC{Float64,Int}`: the edge weight matrix in a sparse format
* `xy::Matrix{Float64}`: the coordinates of the vertices (could be >= 3D)
* `f::Matrix{Float64}`: the values of the function/data at the vertices; each column is an input signal vector, i.e., number of rows == number of nodes.
* `name::String`: a title/name for the data
* `plotspecs::Vector{Symbol}`: specifications for plotting:
    (most likely may be dropped later thanks to julia/Plots.jl)
                   - symm: use a symmetric colorscale
                   - gray: use grayscale
                   - gray255: plot a grayscale image with color bounds [0,255]
                   - copper: use a copper color scale
                   - notitle: don't display the title
                   - nocolorbar: don't display a colorbar
                   - stem: use a stem plot
                   - CLim[cmin,cmax]: set the dynamic display range to [cmin,cmax]
                   - size25: set the size of the nodes to 25 (or whatever value is specified)
                   - LineWidth2: set the width of the lines to 2 (or whatever value is specified)
                   - LineColorb: set the color of the lines/graph edges to color 'b' (or whatever color is specified)
                   - red/black/blue: make the nodes red, black, or blue
                   - verbatim{{...}}: execute the command in the brackets
"""

type GraphSig
    W::SparseMatrixCSC{Float64,Int} # edge weight matrix
    xy::Matrix{Float64}               # spatial coordinates (could be >= 3D)
    f::Matrix{Float64}                # the data sequences (signals) on the nodes
    name::String                      # name of the GraphSig object
    plotspecs::Vector{Symbol}         # special instructions for plotting
    # The following three are deriable from W, xy, and f:
    length::Int                     # length of the data sequence
    dim::Int                        # the dimension of the node coordinates
    size::Tuple{Int,Int}          # size of all the data specified

    # An inner constructor here.
    function GraphSig(W::SparseMatrixCSC{Float64,Int};
                      xy::Matrix{Float64} = Matrix{Float64}(0, 0),
                      f::Matrix{Float64} = zeros(Base.size(W, 1), 1),
                      name::String = "",
                      plotspecs::Vector{Symbol} = Vector{Symbol}(0),
                      length::Int = Base.size(W, 1),
                      dim::Int = Base.size(xy, 2),
                      size::Tuple{Int,Int} = (Base.size(W, 1), Base.size(W, 2) + Base.size(xy, 2) + Base.size(f, 2)))
#
        if isempty(W)
            error("Cannot setup a graph without an edge weight matrix W!")
        elseif Base.size(W, 1) != Base.size(f, 1)
            error("Signal length != # Nodes!!")
        elseif Base.size(W, 1) != Base.size(xy, 1) && Base.size(xy, 1) != 0
            warn("# Node coordinates != # Nodes. Hence, use null coordinates.")
            new( W, Matrix{Float64}(0, 0), f, name, plotspecs, Base.size(W, 1),
                 0, (Base.size(W, 1), Base.size(W, 2) + Base.size(f, 2)) )
        else # Finally, the full constructor!
            new( W, xy, f, name, plotspecs, Base.size(W, 1), Base.size(xy, 2),
                 (Base.size(W, 1), Base.size(W, 2) + Base.size(xy, 2) + Base.size(f, 2)) )
        end
    end # of the inner constructor
end # of type GraphSig

#
# Other GraphSig related functions
#
"""
    replace_data!(G::GraphSig, newdata::Matrix{Float64})

Replace the data in a GraphSig object

### Input Arguments
* `G::GraphSig`: the input GraphSig object whose data is to be replaced; various fields are updated after execution of this function
* `newdata::Matrix{Float64}`: the new data (must be a matrix)
"""
function replace_data!(G::GraphSig, newdata::Matrix{Float64})

    (rows, cols) = Base.size(newdata)
    if G.length == rows         # update the relevant fields accordingly
        G.f = newdata
        G.size = (rows, Base.size(G.W, 2) + Base.size(G.xy, 2) + cols)
    else
        error("The length of each signal (", rows, ") does not match the number of nodes (", G.length, ")!")
    # elseif G.length == cols # This is a dangerous assumption; so prohibit it.
    #    G.f = newdata'
    # elseif isscalar(newdata) # This does not work anymore unlike MATLAB
    #    G.f = newdata*ones(G.length, 1)
    end
end # of function replace_data!


"""
    (value, sigma) = snr(G1::GraphSig, G2::GraphSig, stdp::Bool = false)

    Compute the SNR between G1 (original signal) and G2 (noisy signal)

### Input Arguments
* `G1::GraphSig`: Original reference graph signal
* `G2::GraphSig`: Restored or noisy graph signal
* `stdp::Bool`: A flag to specify to compute the standard deviation of the difference

### Output Arguments
* `value`: Signal/Noise Ratio in the dB unit
* `sigma`: the standard deviation of the noise (if stdp == true)
"""
function snr(G1::GraphSig, G2::GraphSig, stdp=false)
# In julia, the Frobenius norm of a matrix can be computed via `vecnorm`.
  tmp = vecnorm(G2.f - G1.f)
  if tmp < 10 * eps()
    warn("Two signals are the same, hence SNR = Inf")
    value = Inf
  else
    value = 20 * log10(vecnorm(G1.f) / tmp)
  end

  if stdp
    return value, std(G2.f - G1.f)
  else
    return value
  end

end # of snr


"""
    G = gpath(N, f, name)

Generate a GraphSig object for a 1-D path of length `N`

### Input Arguments
* `N::Int`: the length of the path
* `f::Matrix{Float64}`: the signal to be used
* `name::String`: the name for the GraphSig object

### Output Argument
* `G::GraphSig`: the GraphSig ojbect representing the 1-D path of length `N`
"""
function gpath(N::Int,
               f::Matrix{Float64} = ones(N, 1),
               name::String = "Path of length $(N)")

    # the x coordinates
    xy = Matrix{Float64}(reshape(1:N, N, 1))
    # another option: xy = (1.0:Float64(N))[:,:]

    # the weight matrix
    ev = ones(N - 1)               # subdiagonal entries
    W = spdiagm((ev, ev), (-1, 1)) # size(W) is (N, N)
    # MATLAB:
    # ev = ones(N,1);
    # W = spdiags([ev 0*ev ev], [-1 0 1], N, N);

    # create and return the GraphSig object
    return GraphSig(W, xy = xy, f = f, name = name)

end # of function gpath

"""
    function ggrid(Nx::Int, Ny::Int, connect::Symbol)

Generate a GraphSig object for a 2-D grid that is `Nx` by `Ny`.

### Input Arguments
* `Nx::Int`: the number of points in the x direction
* `Ny::Int`: the number of points in the y direction
* `connect::Symbol`: specifies the connectivity mode of the 2-D grid:
    choices are: :c4 (default), :c8, :full

### Output Argument
* `G::GraphSig`: the GraphSig ojbect representing the 2-D `Nx` by `Ny` grid.
"""
function ggrid(Nx::Int, Ny::Int, connect::Symbol = :c4)

    # the total number of nodes
    N = Nx * Ny

    # the xy coordinates
    # MATLAB:
    # xy = [(1:N)', repmat((1:Ny)',Nx,1)];
    # xy(:,1) = ceil(xy(:,1)/Ny);
    xy =hcat( ceil((1:N) / Ny), repmat(1:Ny, Nx) )

    if connect == :c4         # this is the default, i.e., l; r; u; d.
        # make 1-D weight matrices
        ex = ones(Nx - 1)
        Wx = spdiagm((ex, ex), (-1, 1))
        ey = ones(Ny - 1)
        Wy = spdiagm((ey, ey), (-1, 1))
        # form the 2-D weight matrix
        W = kron(speye(Nx), Wy) + kron(Wx, speye(Ny)) # this keeps a sparse form
    else                # now for 8-connected or fully-connected cases
        W = sparse(pairwise(Euclidean(), xy')) # This is from Distances package.
        W[1:(N + 1):(N * N)] = 10 # set each diagonal entry to 10 to prevent blowup.
        W = W.^-1         # convert dist to affinity = inv. dist.
        if connect == :c8 # 8-connected case incl. diagonal connections
            W[find(W .< 0.7)] = 0
        elseif connect == :full # a complete graph case
            W[1:(N + 1):(N * N)] = 0 # set the diagonal to zero (i.e., no loops)
        else
            error("Cannot understand this connectivity: ", connect)
        end
    end

    return GraphSig(W, xy = xy, f = Matrix{Float64}((1:N)[:, :]), name = "$(Nx) by $(Ny) grid with connectivity mode $(connect)")

end # of function ggrid





"""
    function Adj2InvEuc(G::GraphSig)

Given a GraphSig object 'G' with a binary weight (adjacency) matrix, generate a GraphSig object 'GInvEuc' with an inverse Euclidean distance matrix

### Input Arguments
* `G::GraphSig`: a GraphSig object

### Output Argument
* `G::GraphSig`: a GraphSig object with inverse Euclidean weights
"""
function Adj2InvEuc(G::GraphSig)
    W = deepcopy(G.W)
    xy = deepcopy(G.xy)
    f = deepcopy(G.f)

    #find the edges
    (rows,cols) = Base.findnz(W)

    #remove duplicate points
    for j = length(rows):-1:1
      if j<=length(rows) && norm(xy[rows[j],:]-xy[cols[j],:],2) < 10^3*eps()
        W = W[setdiff(1:end, rows[j]), :]
        W = W[:, setdiff(1:end, rows[j])]
        xy = xy[setdiff(1:end, rows[j]), :]
        f = f[setdiff(1:end, rows[j]),:]
        (rows,cols) = Base.findnz(W)
      end
    end

    for j = 1:length(rows)
      W[rows[j],cols[j]] = 1/norm(xy[rows[j],:]-xy[cols[j],:],2)
    end

    return GraphSig(W, xy = xy, f = f, name = string(G.name, " (inverse Euclidean weight matrix)"), plotspecs = G.plotspecs)
end

"""
    function AddNoise(G::GraphSig; SNR::Float64=Float64(1.24), noisetype::String = "gaussian")

Add noise to the data of a GraphSig object

### Input Arguments
* `G::GraphSig`: a GraphSig object
* `SNR`: the SNR that the noisy signal should have
* `noisetype`: the type of noise: Gaussian (default) or Poisson

### Output Argument
* `G::GraphSig`: the GraphSig object with added noise
* `sigma`: the standard deviation of the noise
"""


function AddNoise(G::GraphSig; SNR::Float64=Float64(1.24), noisetype::String = "gaussian")
  f = deepcopy(G.f)

  if noisetype == "gaussian"
    # generate Gaussian noise
    noise = randn(size(f))

    # scale the noise to the desired SNR level
    sigma = vecnorm(f,2)/10.0^(0.05*SNR)/vecnorm(noise,2)
    noise = sigma*noise

    # generate the noisy signal
    f = f + noise

  elseif noisetype == "poisson"
    # generate Poisson noise
    (rows,cols) = size(f)
    noise = zeros(rows,cols)
    for i = 1:rows
      for j = 1:cols
        noise[i,j] = rand(Poisson(f[i,j])) - f[i,j]
      end
    end

    # scale the noise to the desired SNR level
    sigma = vecnorm(f,2)/10.0^(0.05*SNR)/vecnorm(noise,2)
    noise = sigma*noise

    # generate the noisy signal
    f = f + noise

  end

  return return GraphSig(G.W, xy = G.xy, f = f, name = string("SNR = ",SNR," ", G.name), plotspecs = G.plotspecs)
end

end # of module GraphSignal
