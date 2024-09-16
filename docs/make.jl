#as long as source directory is not accessible through Julia's LOAD_PATH
push!(LOAD_PATH,"../src/")

using Documenter, DenseGillespieAlgorithm

makedocs(sitename="DenseGillespieAlgorithm.jl Documentation",
         pages = [
            "Index" => "index.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/roccminton/DenseGillespieAlgorithm.jl.git",
    devbranch = "main"
)
