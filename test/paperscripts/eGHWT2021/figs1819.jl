# Script to generate Figures 18 and 19
using FileIO, Images, ImageFiltering, JLD2, MultiscaleGraphSignalTransforms, MultivariateStats, Plots, StatsBase
include("auxilaries.jl")

# Set up the resolution and display size
gr(dpi=200, size=(800,600))

# Load the original 512x512 composite texture image and its mastge
textures = FileIO.load("5block.png")
textures = Matrix{Float64}(textures) # convert it to a regular matrix

# Compute the mask image using the Gabor features + PCA
#
# Step 1: Parameter setups
# Orientations
ndir = 2
Θ = [ π/3; 5*π/6 ]

# Spatial frequencies
Ξ = [0.2; 0.3; 0.4; 0.5]
# Convert them to wavelengths
Λ = 1.0 ./ Ξ

# Gaussian bandwidth for Gabor filters
bw = 1.0
σ = 1.0/π * sqrt(log(2)/2) * (2^bw+1)/(2^bw-1)

# Spatial aspect ratio: σx = σ; σy = σ/γ
γ = 1.0

# Step 2: Apply Gabor filters; the same size as the input image.
m, n = size(textures)

K = length(Θ) *length(Λ) # K should be square number
F = zeros(m, n, K) # This is the space for the Gabor filtered images
Fg = zeros(m, n, K) # This is the space for the Gaussian smoothed version of F
kr = zeros(m+1, n+1, K)
ki = zeros(m+1, n+1, K)

k = 1
kl = (0, 0)
for θ in Θ
    for λ in Λ
        kernel = Kernel.gabor(m, n, λ*σ, θ, λ, γ, 0.0)
        kernel[1] ./= 2*π*σ*σ/γ
        kernel[2] ./= 2*π*σ*σ/γ
        global kl = size(kernel[1])
        kr[1:kl[1],1:kl[2],k] = kernel[1]
        ki[1:kl[1],1:kl[2],k] = kernel[2]
        cker = (centered(kernel[1]), centered(kernel[2]))
        F[:,:,k] = sqrt.(imfilter(textures, reflect(cker[1])).^2 + 
            imfilter(textures, reflect(cker[2])).^2)
        # "reflect"ing the kernel is necessary for convolution;
        # otherwise it does the correlation.
        Fg[:,:,k] = imfilter(F[:,:,k], Kernel.gaussian(3.0*λ))
        global k += 1
    end
end 

# Display those Gabor features
plot(K, layout=(length(Θ),length(Λ)), framestyle=:none, legend=false)
for k = 1:K
    heatmap!(kr[226:286,226:286,k], subplot=k, ratio=1, yaxis=:flip, showaxis=false, 
        ticks=false, c=:grays, clims=(minimum(kr), maximum(kr)), colorbar=false)
end
display(current())

# Displaying the imaginary part of Gabor kernels may not be necessary...
#plot(K, layout=(Int(sqrt(K)),Int(sqrt(K))), framestyle=:none, legend=false)
#for k = 1:K
#    heatmap!(ki[:,:,k], subplot=k, ratio=1, yaxis=:flip, showaxis=false, 
#      ticks=false, c=:grays, clims=(minimum(ki), maximum(ki)), colorbar=false)
#end
#display(current())

plot(K, layout=(length(Θ),length(Λ)), framestyle=:none, legend=false)
for k = 1:K
    heatmap!(Fg[:,:,k], subplot=k, ratio=1, yaxis=:flip, showaxis=false, 
        ticks=false, c=:grays, clims=(0, maximum(F)), colorbar=false)
end
display(current())

# Step 3.  Compute the PCA and extract the first principal component
Xtr=reshape(Fg, (m*n, K))
Xtr=Xtr'
dt=StatsBase.fit(ZScoreTransform, Xtr) # Each variable is standardized to have 
# mean 0, stddev 1
Xtr=StatsBase.transform(dt, Xtr)
M = MultivariateStats.fit(PCA, Xtr)
Ytr=MultivariateStats.transform(M, Xtr)
Ytr=Ytr'
Ytr=reshape(Ytr, (m, n, size(Ytr, 2)))
plot(K, layout=(length(Θ),length(Λ)), framestyle=:none, legend=false)
for k = 1:outdim(M)
    heatmap!(Ytr[:,:,k], subplot=k, ratio=1, yaxis=:flip, showaxis=false, 
             ticks=false, c=:grays, clims=(minimum(Ytr), maximum(Ytr)), colorbar=false)
end
display(current())

mask = Ytr[:,:,1]
display(heatmap(mask, ratio=1, c=:grays, yaxis=:flip))
# End of 512 x 512 mask generation.

#
# Subsample both the original and mask images to 128 x 128
#
textures = deepcopy(textures[1:4:end,1:4:end])
mask = deepcopy(mask[1:4:end,1:4:end])
m, n = size(textures) # recapture the subsampled matrix size

#
# Normalize the mask
#
mask = (mask.-minimum(mask))./(maximum(mask)-minimum(mask));

#
# Generate Fig. 18a
#
display(heatmap(textures, ratio=1, yaxis =:flip, showaxis = false, ticks = false,
        color = :grays, colorbar = false))
savefig("textures_orig.pdf")
savefig("textures_orig.png")

#
# Generate Fig. 18b
#
display(heatmap(mask, ratio=1, yaxis =:flip, showaxis = false, ticks = false,
        color = :grays, colorbar = false))
savefig("textures_mask.pdf")
savefig("textures_mask.png")

#
# Now, let's compute the weight matrix based on the Gaussian affinity
# of the small windows (or radius r) around pixels and the pixel locations.
#
# Set up the key parameters
r = 3 # specify radius of neighbors; so far the best
σ = 0.0005 # the best so far with r=3 or r=5

# Do the weight matrix computation
W = image_Gaussian_affinity(mask, r, σ)    

#
# Generate G (GraphSig object) and GP (GraphPart object) using the computed
# weight matrix W, and the subsampled texture image as a graph signal.
# Note that GraphSig requires a matrix data even if it is just a one vector, i.e., f = textures[:] does not work!
G = GraphSig(W, f = reshape(textures, (length(textures), 1)))
GP = partition_tree_fiedler(G)
dmatrix = ghwt_analysis!(G, GP=GP)

#
# Construct or search the specific basis
#
# Haar
BS_haar = bs_haar(GP)
dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)

# Walsh
BS_walsh = bs_walsh(GP)
dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)

# GHWT_c2f
dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)

# GHWT_f2c
dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)

# eGHWT
dvec_eghwt, BS_eghwt = eghwt_bestbasis(dmatrix, GP)

#
# Generate Figure 19a (Approximation error plot)
#
DVEC = [ dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:] ]
T = [ "Graph Haar", "Graph Walsh", "GHWT_c2f ", "GHWT_f2c", "eGHWT" ]
L = [ (:dashdot,:orange), (:dashdot,:blue), (:solid, :red), (:solid, :green), 
      (:solid, :black) ]
approx_error2(DVEC, T, L, 0.5)
savefig("textures_approx_error_r3_sigma0005.pdf")
savefig("textures_approx_error_r3_sigma0005.png")

#
# Generate Figure 19b (Top 9 eGHWT basis vectors)
#
top_vectors_plot3(dvec_eghwt, m, n, BS_eghwt, GP)
savefig("textures_eghwt09_r3_sigma0005.pdf")
savefig("textures_eghwt09_r3_sigma0005.png")
