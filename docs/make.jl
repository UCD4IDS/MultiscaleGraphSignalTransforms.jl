using Documenter, MultiscaleGraphSignalTransforms

makedocs(
    sitename="MultiscaleGraphSignalTransforms.jl",
    format = Documenter.HTML(),
    authors = "Jeff Irion, Haotian Li, Naoki Saito, Yiqun Shao",
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "Functions" => "functions.md",
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/haotian127/MultiscaleGraphSignalTransforms.jl.git"
)
