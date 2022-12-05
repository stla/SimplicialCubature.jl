import Pkg
Pkg.add("Documenter")

push!(LOAD_PATH, "../src/")

using Documenter, SimplicialCubature

makedocs(sitename = "SimplicialCubature.jl")

deploydocs(
    repo = "github.com/stla/SimplicialCubature.jl.git"
)
