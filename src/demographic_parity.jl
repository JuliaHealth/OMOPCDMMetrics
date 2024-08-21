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
        _, true_subjects, _ =  _overlapped_subjects(cohorts_df)

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

export demographic_parity
