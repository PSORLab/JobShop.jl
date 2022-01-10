
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function terminated(j::JobShopProblem)
    if j.status.current_iteration > j.parameter.iteration_limit
        println("Terminated as current_iteration $(j.status.current_iteration) > iteration limit $(j.parameter.iteration_limit).")
        return true
    elseif current_abs_gap(j) <= j.parameter.absolute_tolerance
        println("Terminated as absolute gap $(current_abs_gap(j)) <= absolute tolerance $(j.parameter.absolute_tolerance).")
        return true
    elseif current_rel_gap(j) <= j.parameter.relative_tolerance
        println("Terminated as relative gap $(current_rel_gap(j)) <= relative tolerance $(j.parameter.relative_tolerance).")
        return true
    end
    return false
end

function sequential_solve!(j::JobShopProblem)
    initialize!(j)
    k = 1
    while !terminated(j)
        if solve_subproblem(j, j.Ii[k])
            if use_problem(FeasibilityProblem(), j)
                solve_problem(FeasibilityProblem(), j)
            end
            j.status.current_M += 1
            j.lambd[j.status.current_M] .= j.mult
            if use_problem(StepsizeProblem(), j)
                solve_problem(StepsizeProblem(), j)           
            end
            j.status.prior_norm = j.status.current_norm
            j.status.prior_step = j.status.current_step
            display_iteration(j, k)
            k = mod(current_iteration(j), length(j.Ii)) + 1
            j.status.current_iteration += 1
        end
    end
    return nothing
end
