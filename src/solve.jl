
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

function sequential_solve(d::JobShopProblem)
    for i in d.I
        m, o, g, s = create_subproblem(m, i, λ)
        jsprob.m[i] = m
        jsprob.o[i] = o
        jsprob.g[i] = g
        jsprob.s[i] = s
    end
    # TODO: initialize duals
    while !check_termination!(d)
        j = mod(k, length(d.I)) + 1 
        update_subproblem!(d, j, λ)
        solve_subproblem!(d, j)
        update_search_direction!(d)
    end
    solve_feasibility_problem(d)
end

function parallel_solve(jsprob::JobShopProblem)
end
