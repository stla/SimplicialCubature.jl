import Pkg
Pkg.add("Documenter")

push!(LOAD_PATH, "../src/")

using Documenter, SimplicialCubature

makedocs(sitename = "SimplicialCubature.jl")
