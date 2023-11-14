push!(LOAD_PATH, "../src/")

using Documenter, SimplicialCubature

makedocs(sitename = "SimplicialCubature.jl", modules = [SimplicialCubature])

deploydocs(
    repo = "github.com/stla/SimplicialCubature.jl.git"
)
