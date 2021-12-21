
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function check_termination!(d::JobShopProblem)
    (d.status.current_iteration >= d.parameter.iteration_limit) && return true
    (current_abs_tol(d) >= d.parameter.absolute_tolerance)      && return true
    (current_rel_tol(d) >= d.parameter.relative_tolerance)      && return true
    d.status.current_iteration += 1
    return false
end

function update_norm!(d::JobShopProblem, j)
    for m ∈ d.M, t ∈ d.T
        d.status.current_norm += max(d.s[j][m,t], 0.0)^2
    end
    return nothing
end

function save_solution!(::Subproblem, d::JobShopProblem, j)
    d.status.current_lower_bound = objective(d.m[j])
    return nothing
end

function update_step!(d::JobShopProblem, j)
    @unpack current_lower_bound, current_estimate, current_norm = d.status
    α = d.parameter.alpha_step
    d.status.current_step = α/10.0*(current_estimate - current_lower_bound)/current_norm
    return nothing
end

function sequential_solve(d::JobShopProblem)
    @unpack M, T = d
    d.status.current_estimate = d.parameter.start_estimate
    d.status.prior_norm = d.status.current_norm = d.parameter.start_norm
    d.status.prior_step = d.status.current_step = d.parameter.start_step
    for i in d.I
        m, o, g, s = create_problem(Subproblem(), m, i, λ)
        jsprob.m[i] = m
        jsprob.o[i] = o
        jsprob.g[i] = g
        jsprob.s[i] = s
    end
    while !check_termination!(d)
        j = mod(k, length(d.I)) + 1 
        update_subproblem!(Subproblem(), d, j, λ)
        optimize!(d.m[j])
        d.status.solve_time += solve_time(d.m[j])
        if check_termination(Subproblem(), d.m[j])
            save_solution!(Subproblem(), d, j)
            update_norm!(d, j)
            update_step!(d, j)
            update_search_direction!(d)
            if check_feasible_lambda(d)
                mstep = create_problem(StepsizeProblem(), d)
            end
            d.status.current_iteration += 1
            display_iteration()
        end
    end
    solve_feasibility_problem(d)
end

