module Fairness

    using DataFrames
    using OMOPCDMCohortCreator

    import Combinatorics:
        powerset    
    import Chain:
       @chain 

    function _overlapped_subjects(cohorts, conn)

        required_overlapping_phenotypes = ceil(length(cohorts) / 2)
        intersecting_phenotypes = filter(x -> length(x) == required_overlapping_phenotypes, collect(powerset(cohorts)))

        intersecting_subjects = []
        for overlap in map(overlap_subset -> GetCohortSubjects.(overlap_subset, conn), intersecting_phenotypes)
            subjects = []
            for cohort in overlap
                push!(subjects, cohort.subject_id)
            end
            push!(intersecting_subjects, intersect(subjects...))
        end

        # This is the Y from the paper!
        full_subjects = union(intersecting_subjects...)

    end

    function demographic_parity(cohorts, classes, conn; silver_standard = false)

        if silver_standard == true
            dp = Dict(k => Dict() for k in vcat(cohorts, :silver_standard))
        else
            dp = Dict(k => Dict() for k in cohorts)
        end

        patients = GetDatabasePersonIDs(conn)

        for class in classes
            for cohort in cohorts
                cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

                cohort_df = @chain class(cohort_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :cohort_count)
                end

                total_df = @chain class(patients, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :total_count)
                end

                join_feature = intersect(names(cohort_df), names(total_df))

                df = outerjoin(cohort_df, total_df, on = join_feature, makeunique = true, matchmissing = :equal)

                df.dp = df.cohort_count ./ df.total_count

                push!(dp[cohort], Symbol(join_feature[1]) => df)
            end
        end

        if silver_standard == true
        
            Y = _overlapped_subjects(cohorts, conn)

            study_subjects = filter!(
                x -> in(x, Y), 
                GetCohortSubjects(cohorts, conn).subject_id
            ) |> unique

            for class in classes

                study_df = @chain class(study_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :study_count)
                end

                total_df = @chain class(patients, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :total_count)
                end

                join_feature = intersect(names(study_df), names(total_df))

                df = outerjoin(study_df, total_df, on = join_feature, makeunique = true, matchmissing = :equal)

                df.dp = df.study_count ./ df.total_count

                push!(dp[:silver_standard], Symbol(join_feature[1]) => df)
            end
        end

        return dp 

    end


    #=
    #
    # What is now clear to me is that demographic parity is really a more generalized proportion.
    # In a similar manner, this can be turned into a function that is more of a "class proportion"
    # function. Explore ways to restrict just to cohort or to overall database. Additionally, 
    # notice how this could be used to calculate prevalence very rapidly as well as demographic
    # parity.
    #
    =#


    function equality_of_opportunity(cohorts, classes, conn; silver_standard = false)

        if silver_standard == true
            eoo = Dict(k => Dict() for k in vcat(cohorts, :silver_standard))
        else
            eoo = Dict(k => Dict() for k in cohorts)
        end

        Y = _overlapped_subjects(cohorts, conn)
        study_subjects = filter!(
            x -> in(x, Y), 
            GetCohortSubjects(cohorts, conn).subject_id
        ) |> unique

        for cohort in cohorts

            cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

            for class in classes 

                cohort_df = @chain class(cohort_subjects, conn) begin
                    filter(row -> in(row[1], Y), _)
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :cohort_count)
                end

                study_df = @chain class(study_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :study_count)
                end

                join_feature = intersect(names(cohort_df), names(study_df))

                df = outerjoin(cohort_df, study_df, on = join_feature, makeunique = true, matchmissing = :equal)

                df.eoo = df.cohort_count ./ df.study_count

                push!(eoo[cohort], Symbol(join_feature[1]) => df)

            end

        end

        if silver_standard == true
        
            Y = _overlapped_subjects(cohorts, conn)

            study_subjects = filter!(
                x -> in(x, Y), 
                GetCohortSubjects(cohorts, conn).subject_id
            ) |> unique

            for class in classes

                df = @chain class(study_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :study_count)
                end

                df.eoo = df.study_count ./ df.study_count

                push!(eoo[:silver_standard], Symbol(names(df)[1]) => df)
            end
        end

        return eoo
        
    end

    function predictive_rate_parity(cohorts, classes, conn; silver_standard = false)

        if silver_standard == true
            prp = Dict(k => Dict() for k in vcat(cohorts, :silver_standard))
        else
            prp = Dict(k => Dict() for k in cohorts)
        end

        Y = _overlapped_subjects(cohorts, conn)

        for cohort in cohorts

            cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

            for class in classes 

                overlapped_cohort_df = @chain class(cohort_subjects, conn) begin
                    filter(row -> in(row[1], Y), _)
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :overlapped_count)
                end

                cohort_df = @chain class(cohort_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :cohort_count)
                end

                join_feature = names(overlapped_cohort_df)[1]

                df = outerjoin(overlapped_cohort_df, cohort_df, on = join_feature, makeunique = true, matchmissing = :equal)

                df.prp = df.overlapped_count ./ df.cohort_count

                push!(prp[cohort], Symbol(join_feature) => df)

            end

        end

        if silver_standard == true
        
            Y = _overlapped_subjects(cohorts, conn)

            study_subjects = filter!(
                x -> in(x, Y), 
                GetCohortSubjects(cohorts, conn).subject_id
            ) |> unique

            for class in classes

                df = @chain class(study_subjects, conn) begin
                    #=

                    I don't like this 2:end business but it is the
                    closest generalized way I can think of to drop
                    the central axis of a cohort (i.e. the derived
                    subject ID from a given phen definition).
                    The reason for that is because the actual 
                    meaning of a subject ID can vary depending on
                    what cohort has been designed.

                    =#
                    groupby(_, names(_)[2:end])
                    combine(_, nrow => :study_count)
                end

                df.prp = df.study_count ./ df.study_count

                push!(prp[:silver_standard], Symbol(names(df)[1]) => df)
            end
        end

        return prp
        
    end

end
