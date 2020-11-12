using Documenter, Coluna

makedocs(
    modules = [Coluna],
    checkdocs = :exports,
    sitename = "Coluna User Guide",
    format = Documenter.HTML(),
    strict = true,
    pages = Any[
        "Introduction"   => "index.md",
        "Manual" => Any[
            "Getting started"   => "user/start.md",
            "Callbacks"   => "user/callbacks.md"
        ],
        "Reference" => Any[
            "Algorithms" => "dev/algorithms.md",
            "Formulation" => "dev/formulation.md",
            "Reformulation" => "dev/reformulation.md",
            "TODO" => "dev/todo.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/atoptima/Coluna.jl.git",
    push_preview = true
)
