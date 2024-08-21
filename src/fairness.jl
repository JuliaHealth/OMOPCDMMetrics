module Fairness

    using DataFrames
    using OMOPCDMCohortCreator
    import Base:
        Fix2

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

    """
    function demographic_parity(
        cohorts, 
        funcs, 
        conn; 
        labels = false, 
        silver = false, 
        reference_subjects = "", 
        subset_length = 10000
    )


    """
    function demographic_parity(
        cohorts, 
        funcs, 
        conn; 
        labels = false, 
        silver = false, 
        reference_subjects = "", 
        subset_length = 10000
    )
        if labels == true
            _demographic_parity(cohorts, funcs, conn,
            reference_subjects, 
            subset_length,
            silver)
        else
            _demographic_parity(cohorts, funcs, conn, reference_subjects, 
            subset_length)
        end
    end

    """
    _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, subset_length, silver)

    """
    function _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, subset_length, silver)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        cohorts_df = GetCohortSubjects(cohorts, conn)

        subsets = _subset_subjects(reference_subjects, subset_length)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, Not(:count_denom)) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        num = DataFrame()
        for cohort_idx in unique(cohorts_df.cohort_definition_id)
            subjects = filter(row -> row.cohort_definition_id == cohort_idx, cohorts_df).subject_id
            subsets = _subset_subjects(subjects, subset_length)
            for sub in subsets
                vals = _counter_reducer(sub, :count_num, _funcs)
                vals.cohort_definition_id .= cohort_idx
                num = vcat(num, vals)
            end
        end

        if silver == true
            _, true_subjects, _ =  _overlapped_subjects(cohorts, conn)

            subsets = _subset_subjects(true_subjects, subset_length)

            silver = DataFrame()
            for sub in subsets
                silver = vcat(silver, _counter_reducer(sub, :count_num, _funcs))
            end

            silver.cohort_definition_id .= :silver
            silver = groupby(silver, Not(:count_num)) |> 
            x -> combine(x, :count_num => sum => :count_num)

            num = vcat(num, silver)

        end

        num = groupby(num, Not(:count_num)) |> 
        x -> combine(x, :count_num => sum => :count_num)

        dps = outerjoin(num, denom; on = names(num)[1:end-2] .|> 
        x -> Symbol(x) => Symbol(x)) |>
        x -> coalesce.(x, 0)

        dps.demographic_parity = dps.count_num ./ dps.count_denom

        return dps
    end

    """
    _demographic_parity(cohorts::DataFrame, funcs, conn, reference_subjects, subset_length, silver)

    """
    function _demographic_parity(cohorts::DataFrame, funcs, conn, reference_subjects, subset_length, silver)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        println("Getting ready to run")
        println(_funcs[1](1))

        cohorts_df = cohorts
        cohorts = cohorts.subject_id

        subsets = _subset_subjects(reference_subjects, subset_length)

        denom = DataFrame()
        for (idx, sub) in enumerate(subsets)
            println("Denom Iteration is running: $idx")
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, Not(:count_denom)) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        println("Calculating Numerator")
        num = DataFrame()
        for cohort_idx in unique(cohorts_df.cohort_definition_id)
            subjects = filter(row -> row.cohort_definition_id == cohort_idx, cohorts_df).subject_id
            subsets = _subset_subjects(subjects, subset_length)
            for sub in subsets
                vals = _counter_reducer(sub, :count_num, _funcs)
                vals.cohort_definition_id .= cohort_idx
                num = vcat(num, vals)
            end
        end

        println("Skipping Silver Calculation")
        if silver == true
            _, true_subjects, _ =  _overlapped_subjects(cohorts, conn)

            subsets = _subset_subjects(true_subjects, subset_length)

            silver = DataFrame()
            for sub in subsets
                silver = vcat(silver, _counter_reducer(sub, :count_num, _funcs))
            end

            silver.cohort_definition_id .= :silver
            silver = groupby(silver, Not(:count_num)) |> 
            x -> combine(x, :count_num => sum => :count_num)

            num = vcat(num, silver)

        end

        num = groupby(num, Not(:count_num)) |> 
        x -> combine(x, :count_num => sum => :count_num)

        dps = leftjoin(num, denom; 
            on = names(num)[1:end-2] .|> 
            x -> Symbol(x) => Symbol(x)
        )

        dps.demographic_parity = dps.count_num ./ dps.count_denom

        return dps
    end

    """
    _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, subset_length)

    """
    function _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, subset_length)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        cohorts = GetCohortSubjects(cohorts, conn).subject_id

        subsets = _subset_subjects(reference_subjects, subset_length)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, names(denom)[1:end-1]) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        subsets = _subset_subjects(cohorts, subset_length)

        num = DataFrame()
        for sub in subsets
            num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
        end

        num = groupby(num, names(num)[1:end-1]) |> 
        x -> combine(x, :count_num => sum => :count_num)

        dps = outerjoin(num, denom; on = names(num)[1:end-1] .|> 
        x -> Symbol(x) => Symbol(x)) |>
        x -> coalesce.(x, 0)

        dps.demographic_parity = dps.count_num ./ dps.count_denom

        return dps
    end

    """
    _demographic_parity(cohorts::DataFrame, funcs, conn, reference_subjects, subset_length)


    """
    function _demographic_parity(cohorts::DataFrame, funcs, conn, reference_subjects, subset_length)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        println("Getting ready to run")
        println(_funcs[1](1))

        cohorts = cohorts.subject_id

        subsets = _subset_subjects(reference_subjects, subset_length)

        denom = DataFrame()
        for (idx, sub) in enumerate(subsets)
            @info "Denom Iteration $idx"
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, names(denom)[1:end-1]) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        subsets = _subset_subjects(cohorts, subset_length)

        num = DataFrame()
        for (idx, sub) in enumerate(subsets)
            @info "Num Iteration $idx"
            num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
        end

        num = groupby(num, names(num)[1:end-1]) |> 
        x -> combine(x, :count_num => sum => :count_num)

        dps = outerjoin(num, denom; on = names(num)[1:end-1] .|> 
        x -> Symbol(x) => Symbol(x)) |>
        x -> coalesce.(x, 0)

        dps.demographic_parity = dps.count_num ./ dps.count_denom

        return dps
    end

    """
    equality_of_opportunity(cohorts, funcs, conn; reference_subjects = "", subset_length = 10000)

    """
    function equality_of_opportunity(cohorts, funcs, conn; reference_subjects = "", subset_length = 10000)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        study_subjects, true_subjects, false_subjects =  _overlapped_subjects(cohorts, conn)

        subsets = _subset_subjects(true_subjects, subset_length)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, names(denom)[1:end-1]) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        eoo = DataFrame()
        for cohort_idx in cohorts

            cohort = GetCohortSubjects(cohort_idx, conn)
            cohort = filter(row -> in(row.subject_id, true_subjects), cohort)

            subsets = _subset_subjects(cohort.subject_id, subset_length)

            num = DataFrame()
            for sub in subsets
                num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
            end

            num = groupby(num, names(num)[1:end-1]) |> 
            x -> combine(x, :count_num => sum => :count_num)

            cohort = outerjoin(num, denom; on = names(num)[1:end-1] .|> 
            x -> Symbol(x) => Symbol(x)) |>
            x -> coalesce.(x, 0)

            cohort.equality_of_opportunity = cohort.count_num ./ cohort.count_denom

            cohort.cohort_definition_id = ones(Int, nrow(cohort)) .* cohort_idx
            eoo = vcat(eoo, cohort)
        end

        return eoo
    end

    """
    predictive_rate_parity(cohorts, funcs, conn; reference_subjects = "", subset_length = 10000)

    """
    function predictive_rate_parity(cohorts, funcs, conn; reference_subjects = "", subset_length = 10000)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        study_subjects, true_subjects, false_subjects =  _overlapped_subjects(cohorts, conn)

        prp = DataFrame()
        for cohort_idx in cohorts

            cohort = GetCohortSubjects(cohort_idx, conn)
            true_cohort = filter(row -> in(row.subject_id, true_subjects), cohort)
            false_cohort = filter(row -> in(row.subject_id, false_subjects), cohort)
            
            subsets = _subset_subjects(true_cohort.subject_id, subset_length)

            num = DataFrame()
            for sub in subsets
                num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
            end
            
            subsets = _subset_subjects(false_cohort.subject_id, subset_length)

            if !isempty(subsets)
                false_denom = DataFrame()
                for sub in subsets
                    false_denom = vcat(false_denom, _counter_reducer(sub, :count_num, _funcs))
                end
                denom = vcat(num, false_denom)
                denom = groupby(denom, names(denom)[1:end-1]) |> 
                x -> combine(x, :count_num => sum => :count_denom)
            else 
                denom = num
                denom = groupby(denom, names(denom)[1:end-1]) |> 
                x -> combine(x, :count_num => sum => :count_denom)
            end

            num = groupby(num, names(num)[1:end-1]) |> 
            x -> combine(x, :count_num => sum => :count_num)

            cohort = outerjoin(num, denom; on = names(num)[1:end-1] .|> 
            x -> Symbol(x) => Symbol(x)) |>
            x -> coalesce.(x, 0)

            cohort.predictive_rate_parity = cohort.count_num ./ cohort.count_denom

            cohort.cohort_definition_id = ones(Int, nrow(cohort)) .* cohort_idx
            prp = vcat(prp, cohort)
        end

        return prp 
    end

    export demographic_parity, equality_of_opportunity,predictive_rate_parity

end
