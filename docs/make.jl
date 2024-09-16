#as long as source directory is not accessible through Julia's LOAD_PATH
push!(LOAD_PATH,"../src/")

using Documenter, DenseGillespieAlgorithm

makedocs(sitename="My Documentation")
