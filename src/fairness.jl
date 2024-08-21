module Fairness

    using DataFrames
    using OMOPCDMCohortCreator
    import Base:
        Fix2

    include("./helpers.jl")
    include("./demographic_parity.jl")
    include("./equality_of_opportunity.jl")
    include("./predictive_rate_parity.jl")

end
