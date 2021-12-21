
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function terminated(d::JobShopProblem)
    (d.status.current_iteration >= d.parameter.iteration_limit) && return true
    (current_abs_gap(d) >= d.parameter.absolute_tolerance)      && return true
    (current_rel_gap(d) >= d.parameter.relative_tolerance)      && return true
    d.status.current_iteration += 1
    return false
end

function update_norm_step!(d::JobShopProblem, j)
    @unpack current_lower_bound, current_estimate, current_norm = d.status
    for m ∈ d.M, t ∈ d.T
        d.status.current_norm += max(d.s[j][m,t], 0.0)^2
    end
    α = d.parameter.alpha_step
    d.status.current_step = α/10.0*(current_estimate - current_lower_bound)/current_norm
    return nothing
end

function save_solution!(::Subproblem, d::JobShopProblem, j)
    d.status.lower_bound[d.status.current_iteration] = objective(d.m[j])
    d.mult .+= d.status.current_step .* value.(d.s[j])
    d.mult .= max.(d.mult, 0.0)

    d.slack_value .= value.(d.s[j])
    d.vp_value    .= value.(d.vp[j])
    d.sbTimeI1    .= value.(bTimeI1[j])
    d.sdelta1     .= value.(delta1[j])
    d.sbTimeI2    .= value.(m.bTimeI2[j])

    d.sTard1 .= value.(Tard1)
    d.sTard2 .= value.(Tard2)
    d.sbTime1 .= value.(bTime1)
    
    return nothing
end

function sequential_solve(d::JobShopProblem)
    # initialize subproblems & user paramters
    @unpack M, T = d
    d.status.start_time = time()
    d.status.current_estimate = d.parameter.start_estimate
    d.status.prior_norm = d.status.current_norm = d.parameter.start_norm
    d.status.prior_step = d.status.current_step = d.parameter.start_step
    d.status.current_M = d.parameter.start_M
    for i in d.I
        m, o, g, s = create_problem(Subproblem(), d, i, λ)
        jsprob.m[i] = m
        jsprob.o[i] = o
        jsprob.g[i] = g
        jsprob.s[i] = s
    end
    create_problem!(FeasibilityProblem(), d)

    # begin main solution loop
    while !terminated(d)
  
        j = mod(d.status.current_iteration , length(d.I)) + 1 
        d.status.solve_time += update_solve!(Subproblem(), d, j, λ)
  
        if valid_solve(Subproblem(), d.m[j])
            
            save_solution!(Subproblem(), d, j)
            update_norm_step!(d, j)

            if use_problem(FeasibilityProblem(), d)
                d.status.heurestic_time += update_solve!(FeasibilityProblem(), d)
                if valid_solve(FeasibilityProblem(), d.feasibility_model)
                    d.status.upper_bound[k] = objective_value(d.feasibility_model)
                end
            end

            d.status.current_M += 1
            
            if use_problem(StepsizeProblem(), d)
                lambda[M] .= mult
      
                new_maxest = current_step(d)*current_norm(d)/alpha_step(d) + lower_bound(d)
                (d.status.maxest < new_maxest)   && (d.status.maxest = new_maxest)
                (d.status.maxest > d.status.est) && (d.status.maxest = d.status.est)

                d.stepsize_model = create_problem(StepsizeProblem(), d, lambda, M, sstep)
                if (d.status.current_M > 5) && (d.status.current_iteration > 50) && (d.status.est > lower_bound(d))
                    if !valid_solve(StepsizeProblem(), d.stepsize_model) || d.status.current_M >= 50000
                        d.status.est = d.status.maxest 
                        d.status.current_step /= 10
                        d.status.current_M = 1
                        d.status.maxest = -100000
                    end
                end              
            end
            
            d.status.current_norm = d.status.prior_norm
            d.status.current_step = d.status.prior_step
            d.status.current_iteration += 1
            display_iteration(d)
        end
    end
    return
end

