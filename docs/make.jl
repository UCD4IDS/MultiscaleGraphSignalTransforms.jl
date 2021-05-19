using Documenter, MultiscaleGraphSignalTransforms

makedocs(
    sitename="MultiscaleGraphSignalTransforms.jl",
    format = Documenter.HTML(collapselevel = 1),
    authors = "Jeff Irion, Haotian Li, Naoki Saito, Yiqun Shao",
    clean = true,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "1D Path" => "examples/P64.md",
            "Sunflower Graph" => "examples/Sunflower.md",
        ],
        "Functions" => [
            "Recursive Graph Partitioning" => "functions/Partition.md",
            "Hierarchical Graph Laplacian Eigen Transform" => "functions/HGLET.md",
            "Generalized Haar-Walsh Transform" => "functions/GHWT.md",
            "extended Generalized Haar-Walsh Transform" => "functions/eGHWT.md",
            "Natural Graph Wavelet Dictionaries" => "functions/NGWD.md",
            "Utils" => "functions/Utils.md",
        ],
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/haotian127/MultiscaleGraphSignalTransforms.jl.git"
)
