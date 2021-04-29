module pSGWT

using PyCall

export sgwt_frame, sgwt_filter_banks

const pygsp = PyNULL()

function __init__()
    copy!(pygsp, pyimport_conda("pygsp", "pygsp", "conda-forge"))
end

function init_Graph(W)
    G = pygsp.graphs.Graph(PyReverseDims(W))
    return G
end

"""
    sgwt_transform(loc, nf, W)

SGWT\\_TRANSFORM perform the SGWT transform with MexicanHat filter as a python wrapper.
See https://pygsp.readthedocs.io/en/stable/tutorials/wavelet.html

# Input Arguments
- `loc::Int64`: the vertex index that the output wavelet centered at.
- `W::Matrix{Float64}`: weighted adjacency matrix.
- `nf::Int64`: default is 6. Number of filters.

# Output Argument
- `frame::Array{Float64, 3}`: a N x N x nf matrix.

"""
function sgwt_frame(W; nf = 6)
    G = init_Graph(W)
    G.estimate_lmax()
    g = pygsp.filters.MexicanHat(G, Nf = nf)
    # np = pyimport("numpy")
    # s = np.zeros(G.N); s[loc] = 1
    # s = g.filter(s, method="chebyshev")
    return g.compute_frame()
end

function sgwt_filter_banks(W, 𝛌; nf = 6)
    G = init_Graph(W)
    G.estimate_lmax()
    g = pygsp.filters.MexicanHat(G, Nf = nf)
    np = pyimport("numpy")
    e = np.array(𝛌)
    return g.evaluate(e)
end

end
