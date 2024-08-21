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

export predictive_rate_parity
