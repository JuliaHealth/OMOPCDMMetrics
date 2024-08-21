"""
_counter_reducer(sub, count_name, funcs)




"""
function _counter_reducer(sub, count_name, funcs)
    for fun in funcs
        sub = fun(sub)
    end
    sub[:, Not(:person_id)] |>
    x -> groupby(x, names(x)) |> 
    x -> combine(x, nrow => count_name)
end

"""
_subset_subjects(vec::Vector, subset_length::Int; full_subset::Bool = true)

Internal function that accepts a vector, subsets it, and returns a vector of vectors (the subsets of the original vector).

# Arguments:

- `vec::Vector` - a vector of elements to subset

- `subset_length::Int` - how long each subset of the vector should be

# Keywork Arguments:

- `full_subset::Bool` - completely subsets the input vector even if the input vector is not divisible by the given `subset_length` (Default: `true`)

# Returns:

`subsets::Vector{Vector}` - a vector containing subsets of the original vector, `vec`, with each subset of the desired `subset_length`. (NOTE: If `full_subset` is set to `true`, the last element of this vector may be smaller then the `subset_length`). 

"""
function _subset_subjects(vec::Vector, subset_length::Int; full_subset::Bool = true)
    vec_length = size(vec)[1]
    subsets = []
    for i in 1:subset_length:vec_length
        if i + subset_length > vec_length
            if full_subset == true
                push!(subsets, vec[i:end])
            end
        else
            push!(subsets, vec[i:i+subset_length])
        end
    end

    return subsets
end

"""
_overlapped_subjects(cohorts::Vector, conn; patient_overlaps::Symbol = :default)

Internal function which accepts patient cohorts and determines total unique subjects between cohorts, which patients overlap between cohorts, and which patients do not overlap between cohorts.

# Arguments:

- `cohorts::Vector` - a vector of cohort IDs

- `conn` - database connection using DBInterface

# Keyword Arguments:

- `patient_overlaps::Symbol` - determines how subjects are defined as overlapping between cohorts. Available options: `:default`. (NOTE: See References for details on options)

# Returns

A triple in the order, `total_subjects, intersecting_subjects, nonintersecting_subjects` and detailed as follows:

- `total_subjects::Vector` - the total unique subject IDs between all cohorts

- `intersecting_subjects::Vector` - the subject IDs who overlap between cohorts

- `nonintersecting_subjects::Vector` - the subject IDs who do not overlap between cohorts

# References: 

The `:default` algorithm to calculate patient overlaps via `patient_overlaps` is from _T. Y. Sun, S. Bhave, J. Altosaar, and N. Elhadad, “Assessing Phenotype Definitions for Algorithmic Fairness,” arXiv:2203.05174 [cs, q-bio], Mar. 2022, Accessed: Apr. 29, 2022. [Online]. Available: http://arxiv.org/abs/2203.05174_.

"""
function _overlapped_subjects(cohorts, conn; patient_overlaps = :default)

    if patient_overlaps == :default
        required_overlapping_phenotypes = ceil(length(cohorts) / 2)
    end

    subjects = GetCohortSubjects(cohorts, conn)
    subjects.count = [count(==(subject), subjects.subject_id) for subject in subjects.subject_id]

    intersecting_pop = filter(row -> row.count >= required_overlapping_phenotypes, subjects)

    nonintersecting_pop = filter(row -> row.count < required_overlapping_phenotypes, subjects)

    intersecting_subjects = unique(intersecting_pop.subject_id)
    nonintersecting_subjects = unique(nonintersecting_pop.subject_id)
    total_subjects = vcat(intersecting_subjects, nonintersecting_subjects)

    return total_subjects, intersecting_subjects, nonintersecting_subjects 

end
