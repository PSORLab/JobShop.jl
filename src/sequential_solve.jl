
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function terminated(j::JobShopProblem)
    if j.status.current_iteration > j.parameter.iteration_limit
        println("Terminated as current_iteration $(j.status.current_iteration) > iteration limit $(j.parameter.iteration_limit).")
        return true
    elseif j.status.feasible_problem_found
        println("Feasible problem found.")
        return true
    end
    return false
end

function update_stepsize!(jsp::JobShopProblem)
    @unpack prior_norm, current_norm, prior_step, current_iteration = jsp.status
    if current_iteration > 20
        jsp.status.current_step = 0.995*prior_step*sqrt(prior_norm/current_norm)
    end
    return nothing
end

function update_multiplier!(jsp::JobShopProblem)
    @unpack mult, sslackk, MachineType, T = jsp
    @unpack current_step = jsp.status
    for mi in MachineType, t in T
        jsp.mult[mi,t] = max(mult[mi,t] + current_step*sslackk[mi,t], 0.0)
    end
    return nothing
end

function update_penalty!(jsp::JobShopProblem)
    # TODO
    return nothing
end

function sequential_solve!(jsp::JobShopProblem)
    initialize!(jsp)
    k = 0
    while !terminated(jsp)
        if solve_subproblem(jsp, jsp.Ii[k + 1])
            update_stepsize!(jsp)
            update_multiplier!(jsp)
            update_penalty!(jsp)
            if use_problem(FeasibilityProblem(), jsp)
                solve_problem(FeasibilityProblem(), jsp)
            end
            jsp.status.prior_norm = jsp.status.current_norm
            jsp.status.prior_step = jsp.status.current_step
            display_iteration(jsp, k + 1)
            k = mod(k + 1, length(jsp.Ii))
            jsp.status.current_iteration += 1
        end
    end
    return nothing
end
