module Fairness

    import Base:
        Fix2

    import DataStructures:
        counter

    using DataFrames
    using OMOPCDMCohortCreator

    include("./helpers.jl")
    include("./demographic_parity.jl")
    include("./equality_of_opportunity.jl")
    include("./predictive_rate_parity.jl")

end
