module Fairness

    using DataFrames
    using OMOPCDMCohortCreator
    import Base:
        Fix2

    function _counter_reducer(sub, count_name, funcs)
        for fun in funcs
            sub = fun(sub)
        end
        sub[:, Not(:person_id)] |>
        x -> groupby(x, names(x)) |> 
        x -> combine(x, nrow => count_name)
    end

    function _subset_subjects(vec, process_size)
        subsets = []
        for i in 1:process_size:size(vec)[1]
            if i + process_size > size(vec)[1]
                push!(subsets, vec[i:end])
            else
                push!(subsets, vec[i:i+process_size])
            end
        end

        return subsets
    end

    function _overlapped_subjects(cohorts, conn)
    
        # Requirement is based on majority of phenotypes provided
        required_overlapping_phenotypes = ceil(length(cohorts) / 2)

        subjects = GetCohortSubjects(cohorts, conn)
        subjects.count = [count(==(subject), subjects.subject_id) for subject in subjects.subject_id]

        intersecting_pop = filter(row -> row.count >= required_overlapping_phenotypes, subjects)

        nonintersecting_pop = filter(row -> row.count < required_overlapping_phenotypes, subjects)

        intersecting_subjects = unique(intersecting_pop.subject_id)
        nonintersecting_subjects = unique(nonintersecting_pop.subject_id)
        total_subjects = vcat(intersecting_subjects, nonintersecting_subjects)

        return total_subjects, intersecting_subjects, nonintersecting_subjects 

    end

    function demographic_parity(
        cohorts, 
        funcs, 
        conn; 
        labels = false, 
        silver = false, 
        reference_subjects = "", 
        process_size = 10000
    )
        if labels == true
            _demographic_parity(cohorts, funcs, conn,
            reference_subjects, 
            process_size,
            silver)
        else
            _demographic_parity(cohorts, funcs, conn, reference_subjects, 
            process_size)
        end
    end

    function _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, process_size, silver)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        cohorts_df = GetCohortSubjects(cohorts, conn)

        subsets = _subset_subjects(reference_subjects, process_size)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, Not(:count_denom)) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        num = DataFrame()
        for cohort_idx in unique(cohorts_df.cohort_definition_id)
            subjects = filter(row -> row.cohort_definition_id == cohort_idx, cohorts_df).subject_id
            subsets = _subset_subjects(subjects, process_size)
            for sub in subsets
                vals = _counter_reducer(sub, :count_num, _funcs)
                vals.cohort_definition_id .= cohort_idx
                num = vcat(num, vals)
            end
        end

        if silver == true
            _, true_subjects, _ =  _overlapped_subjects(cohorts, conn)

            subsets = _subset_subjects(true_subjects, process_size)

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

    function _demographic_parity(cohorts::Vector{<:Any}, funcs, conn, reference_subjects, process_size)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        cohorts = GetCohortSubjects(cohorts, conn).subject_id

        subsets = _subset_subjects(reference_subjects, process_size)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, names(denom)[1:end-1]) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        subsets = _subset_subjects(cohorts, process_size)

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

    function _demographic_parity(cohorts::DataFrame, funcs, conn, reference_subjects, process_size)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        if isempty(reference_subjects)
            reference_subjects = GetDatabasePersonIDs(conn)
        end

        cohorts = cohorts.subject_id

        subsets = _subset_subjects(reference_subjects, process_size)

        denom = DataFrame()
        for sub in subsets
            denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
        end

        denom = groupby(denom, names(denom)[1:end-1]) |> 
        x -> combine(x, :count_denom => sum => :count_denom)

        subsets = _subset_subjects(cohorts, process_size)

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

    function equality_of_opportunity(cohorts, funcs, conn; reference_subjects = "", process_size = 10000)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        study_subjects, true_subjects, false_subjects =  _overlapped_subjects(cohorts, conn)

        subsets = _subset_subjects(true_subjects, process_size)

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

            subsets = _subset_subjects(cohort.subject_id, process_size)

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

    function predictive_rate_parity(cohorts, funcs, conn; reference_subjects = "", process_size = 10000)

        _funcs = [Fix2(fun, conn) for fun in funcs]

        study_subjects, true_subjects, false_subjects =  _overlapped_subjects(cohorts, conn)

        prp = DataFrame()
        for cohort_idx in cohorts

            cohort = GetCohortSubjects(cohort_idx, conn)
            true_cohort = filter(row -> in(row.subject_id, true_subjects), cohort)
            false_cohort = filter(row -> in(row.subject_id, false_subjects), cohort)
            
            subsets = _subset_subjects(true_cohort.subject_id, process_size)

            num = DataFrame()
            for sub in subsets
                num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
            end
            
            subsets = _subset_subjects(false_cohort.subject_id, process_size)

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
