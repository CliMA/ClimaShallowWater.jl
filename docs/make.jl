using Documenter, DocumenterCitations
using ClimaShallowWater

bib = DocumenterCitations.CitationBibliography(joinpath(@__DIR__, "refs.bib"))

makedocs(
    bib,
    sitename = "ClimaShallowWater",
    format = Documenter.HTML(),
    modules = [ClimaShallowWater],
    pages = [
        "index.md",
        "testcases.md",
        "api.md",
        "references.md",
    ]
)

Documenter.deploydocs(
    repo = "github.com/CliMA/ClimaShallowWater.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
