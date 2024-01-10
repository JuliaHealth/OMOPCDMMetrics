module OMOPCDMCohortMetrics

using DataFrames
using OMOPCDMCohortCreator

import Base.Iterators:
    product

include("./fairness.jl")

function pairwise_overlapping_subjects(cohorts, conn)

    if isa(cohorts, Vector)
        nothing
    else 
        cohorts = [cohorts]
    end

    pairwise_matrix = zeros(length(cohorts), length(cohorts))

    cohort_pairs = product(repeat([cohorts], 2)...) 

    for cohort_pair in cohort_pairs
        groups = groupby(GetCohortSubjects(cohort_pair, conn), :cohort_definition_id)
        if length(groups) == 2
            intersecting_subjects = length(intersect(groups[1].subject_id, groups[2].subject_id))
        else
            intersecting_subjects = length(groups[1].subject_id)
        end
        pairwise_matrix[cohort_pair...] = intersecting_subjects
    end

    return pairwise_matrix

end

end
