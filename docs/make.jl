import Documenter
import LinearAlgebra
import OpenKolmogorovFlow

Documenter.DocMeta.setdocmeta!(
    OpenKolmogorovFlow,
    :DocTestSetup,
    :(using OpenKolmogorovFlow);
    recursive=true,
)

Documenter.makedocs(;
    modules=[OpenKolmogorovFlow],
    sitename="OpenKolmogorovFlow.jl",
    authors="Davide Lasagna",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Davide-Lasagna-s-Lab.github.io/OpenKolmogorovFlow.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "man/getting-started.md",
            "Assumptions" => "man/conventions.md",
            "Fields" => "man/fields.md",
            "Equations" => "man/equations.md",
            "Symmetries" => "man/symmetries.md",
        ],
        "Reference" => "man/api.md",
    ],
    checkdocs=:exports,
)

Documenter.deploydocs(;
    repo="github.com/Davide-Lasagna-s-Lab/OpenKolmogorovFlow.jl.git",
    devbranch="master",
    push_preview=true,
)
