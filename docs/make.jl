cd(@__DIR__) # go into docs folder
import Pkg;
Pkg.activate(@__DIR__);

push!(LOAD_PATH, "../src/")
using Documenter, Literate, PileResponse

# convert tutorial to markdown
Literate.markdown("./src/tutorial.jl", "./src")

# Markdown files to compile to HTML
# (sidebar and table of contents for the documentation)
pages = [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "API" => "api.md"
]

makedocs(;
    pages,
    sitename="PileResponse.jl documentation",
    repo=Remotes.GitHub("antonyorton", "PileResponse.jl")
)

# WARNING: An error occurs if the remote repo has been imported into
#          the base Julia environment. Best to avoid, or use pkg> develop "https://github ..

# comment out when working locally or not updating docs. Only activate when re-deploying a docs update to GitHub
# deploydocs(
#     repo="https://github.com/antonyorton/PileResponse.jl.git",
# )