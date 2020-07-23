using Documenter
using ControleOptimal



makedocs(
modules = [ControleOptimal],
    sitename = "ControleOptimal.jl",
    strict=true,
    authors = "Saloua Naama, Mohamed El Waghf et Rachid ELMontassir",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
            "Accueil" => "index.md"
    ]
    )

deploydocs(repo = "github.com/mathn7/ControlOptimal.git")
