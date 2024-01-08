using OMOPCDMCohortMetrics
using Documenter

DocMeta.setdocmeta!(OMOPCDMCohortMetrics, :DocTestSetup, :(using OMOPCDMCohortMetrics); recursive=true)

makedocs(;
    modules=[OMOPCDMCohortMetrics],
    authors="TheCedarPrince <jacobszelko@gmail.com> and contributors",
    repo="https://github.com/TheCedarPrince/OMOPCDMCohortMetrics.jl/blob/{commit}{path}#{line}",
    sitename="OMOPCDMCohortMetrics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TheCedarPrince.github.io/OMOPCDMCohortMetrics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/TheCedarPrince/OMOPCDMCohortMetrics.jl",
    devbranch="main",
)
