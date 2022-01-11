
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function terminated(j::JobShopProblem)
    if j.status.current_iteration > j.parameter.iteration_limit
        println("Terminated as current_iteration $(j.status.current_iteration) > iteration limit $(j.parameter.iteration_limit).")
        return true
    end
    return false
end

function sequential_solve!(j::JobShopProblem)
    initialize!(j)
    k = 0
    while !terminated(j)
        if solve_subproblem(j, j.Ii[k+1])
            if use_problem(FeasibilityProblem(), j)
                solve_problem(FeasibilityProblem(), j)
            end
            j.status.current_M += 1
            j.lambd[j.status.current_M] .= j.mult
            if j.status.current_iteration > 25
                new_maxest = j.status.current_step*j.status.current_norm/j.parameter.alpha_step_2 + current_lower_bound(j)
                (j.status.maxest < new_maxest)        && (j.status.maxest = new_maxest)
                (j.status.maxest > j.status.estimate) && (j.status.maxest = j.status.estimate)
            end
            if use_problem(StepsizeProblem(), j)
                solve_problem(StepsizeProblem(), j)
            end
            j.status.prior_norm = j.status.current_norm
            j.status.prior_step = j.status.current_step
            display_iteration(j, k+1)
            k += 1
            k = mod(k, length(j.Ii))
            j.status.current_iteration += 1
        end
    end
    return nothing
end
