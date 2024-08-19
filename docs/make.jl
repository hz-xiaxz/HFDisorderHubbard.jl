using HFDisorderHubbard
using Documenter

DocMeta.setdocmeta!(
    HFDisorderHubbard,
    :DocTestSetup,
    :(using HFDisorderHubbard);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
    modules = [HFDisorderHubbard],
    authors = "hz-xiaxz <1806656034@qq.com> and contributors",
    repo = "https://github.com/hz-xiaxz/HFDisorderHubbard.jl/blob/{commit}{path}#{line}",
    sitename = "HFDisorderHubbard.jl",
    format = Documenter.HTML(;
        canonical = "https://hz-xiaxz.github.io/HFDisorderHubbard.jl",
    ),
    pages = [
        "index.md"
        [
            file for file in readdir(joinpath(@__DIR__, "src")) if
            file != "index.md" && splitext(file)[2] == ".md"
        ]
    ],
)

deploydocs(; repo = "github.com/hz-xiaxz/HFDisorderHubbard.jl")
