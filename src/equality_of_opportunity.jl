"""
equality_of_opportunity(cohorts::DataFrame, funcs, conn; reference_subjects = "", subset_length = 10000)

"""
function equality_of_opportunity(cohorts::DataFrame, funcs, conn; reference_subjects = "", subset_length = 10000)

    _funcs = [Fix2(fun, conn) for fun in funcs]

    if isempty(reference_subjects)
        reference_subjects = GetDatabasePersonIDs(conn)
    end

    println("Getting ready to run")
    println(_funcs[1](1))

    study_subjects, true_subjects, false_subjects =  _overlapped_subjects(cohorts)

    subsets = _subset_subjects(true_subjects, subset_length)

    denom = DataFrame()
    for sub in subsets
        denom = vcat(denom, _counter_reducer(sub, :count_denom, _funcs))
    end

    denom = groupby(denom, names(denom)[1:end-1]) |> 
    x -> combine(x, :count_denom => sum => :count_denom)

    eoo = DataFrame()
    for cohort_idx in cohorts.cohort_definition_id |> unique

        cohort = filter(row -> row.cohort_definition_id == cohort_idx, cohorts)
        cohort = filter(row -> in(row.subject_id, true_subjects), cohort)

        subsets = _subset_subjects(cohort.subject_id, subset_length)

        num = DataFrame()
        for sub in subsets
            num = vcat(num, _counter_reducer(sub, :count_num, _funcs))
        end

        num = groupby(num, names(num)[1:end-1]) |> 
        x -> combine(x, :count_num => sum => :count_num)

        cohort = leftjoin(num, denom; 
            on = names(num)[1:end-1] .|> 
            x -> Symbol(x) => Symbol(x)
        )

        cohort.equality_of_opportunity = cohort.count_num ./ cohort.count_denom

        cohort.cohort_definition_id = ones(Int, nrow(cohort)) .* cohort_idx
        eoo = vcat(eoo, cohort)
    end

    return eoo
end

export equality_of_opportunity
