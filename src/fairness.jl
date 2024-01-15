module Fairness

    using DataFrames
    using OMOPCDMCohortCreator

    import Combinatorics:
        powerset    
    import Chain:
       @chain 

    function _overlapped_subjects(cohorts, conn)
    
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

    function demographic_parity(cohorts, classes, conn)

        study, PP, PN = _overlapped_subjects(cohorts, conn)
        
        dps = DataFrame()
        for class in classes
            for cohort in cohorts
                cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

                S = class(study, conn)

                feature_name = names(S)[2]

                for feature in unique(S[:, 2])
                    C = class(cohort_subjects, conn)
                    TP = 
                        filter(row -> row[2] == feature, C) |>
                        filter(row -> in(row[1], PP))

                    FP = 
                        filter(row -> row[2] == feature, C) |>
                        filter(row -> in(row[1], PN))
                
                    N = filter(row -> row[2] == feature, S)

                    dp = (nrow(TP) + nrow(FP)) / nrow(N)

                    push!(dps, Dict(:cohort_definition_id => cohort, Symbol(feature_name) => feature, :dp => dp), cols = :union)

                end
            end
        end

        return dps

    end

    function equality_of_opportunity(cohorts, classes, conn)

        study, PP, PN = _overlapped_subjects(cohorts, conn)
        
        eoos = DataFrame()
        for class in classes
            for cohort in cohorts
                cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

                S = class(study, conn)

                feature_name = names(S)[2]

                for feature in unique(S[:, 2])
                    C = class(cohort_subjects, conn)
                    TP = 
                        filter(row -> row[2] == feature, C) |>
                        filter(row -> in(row[1], PP))

                    P = 
                        filter(row -> row[2] == feature, C)

                    eoo = nrow(TP) / nrow(P)

                    push!(eoos, Dict(:cohort_definition_id => cohort, Symbol(feature_name) => feature, :eoo => eoo), cols = :union)

                end
            end
        end

        return eoos

    end

    function predictive_rate_parity(cohorts, classes, conn)

        study, PP, PN = _overlapped_subjects(cohorts, conn)
        
        prps = DataFrame()
        for class in classes
            for cohort in cohorts
                cohort_subjects = GetCohortSubjects(cohort, conn).subject_id

                S = class(study, conn)

                feature_name = names(S)[2]

                for feature in unique(S[:, 2])
                    C = class(cohort_subjects, conn)
                    TP = 
                        filter(row -> row[2] == feature, C) |>
                        filter(row -> in(row[1], PP))

                    FP = 
                        filter(row -> row[2] == feature, C) |>
                        filter(row -> in(row[1], PN))

                    prp = nrow(TP) / (nrow(TP) + nrow(FP))

                    push!(prps, Dict(:cohort_definition_id => cohort, Symbol(feature_name) => feature, :prp => prp), cols = :union)

                end
            end
        end

        return prps

    end
end
