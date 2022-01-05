
"""
$TYPEDSIGNATURES

Checks termination of the main algorithm.
"""
function terminated(d::JobShopProblem)
    if d.status.current_iteration >= d.parameter.iteration_limit
        println("Terminated as current_iteration $(d.status.current_iteratio) >= iteration limit $(d.parameter.iteration_limit).")
        return true
    elseif current_abs_gap(d) <= d.parameter.absolute_tolerance
        println("Terminated as absolute gap $(current_abs_gap(d)) <= absolute tolerance $(d.parameter.absolute_tolerance).")
        return true
    elseif current_rel_gap(d) <= d.parameter.relative_tolerance
        println("Terminated as relative gap $(current_rel_gap(d)) <= relative tolerance $(d.parameter.relative_tolerance).")
        return true
    end
    return false
end

function update_norm_step!(d::JobShopProblem)
    @unpack alpha_step = d.parameter
    @unpack lower_bound, current_estimate, current_norm = d.status
    @unpack Mi, s = d

    d.status.current_norm += sum(max(sv, 0.0)^2 for sv ∈ s)
    d.status.current_step = alpha_step/10.0*(current_estimate - lower_bound[maximum(keys(lower_bound))])/current_norm
    return nothing
end

function save_solution!(::Subproblem, d::JobShopProblem, m, s, v_p, b1, b2, I)
    @unpack Mi, Tp, J = d
    d.status.lower_bound[d.status.current_iteration] = objective_value(m)
    for m in Mi, t in Tp
        d.s[m,t] = value(s[m,t])
    end
    d.λ .+= d.status.current_step .* d.s
    d.λ .= max.(d.λ, 0.0)
    for i in I, j in J[i]
        d.sb1[i,j] = value(b1[i,j])
    end
    for m in Mi, t in Tp
        d.sv_p[m,t] = value(v_p[m,t])
    end
    #d.bI2  .= value.(bI2)
    return
end

# problems are not preallocated and modified due to RAM considerations 
# (problem 1 - 17 use, ~ 7GB for Pratt Whitney 127 part problem)
"""
$TYPEDSIGNATURES

Runs the main sequential solve algorithm.
"""
function sequential_solve!(d::JobShopProblem)
    # initialize subproblems & user paramters
    @unpack Mi, Tp = d
    d.status.start_time = time()
    d.status.penalty = d.parameter.start_penalty
    d.status.current_estimate = d.parameter.start_estimate
    d.status.prior_norm = d.status.current_norm = d.parameter.start_norm
    d.status.prior_step = d.status.current_step = d.parameter.start_step
    d.status.current_M = d.parameter.start_M
    d.λ = zeros(length(Mi), length(Tp))
    d.status.current_iteration = 0
    d.status.lower_bound[1] = -Inf
    d.status.upper_bound[1] = Inf
    d.stard1 = zeros(length(Mi), length(Tp))
    d.stard2 = zeros(length(Mi), length(Tp), 2)
    d.sslackk = zeros(length(Mi), length(Tp))
    d.sv_p = zeros(length(Mi), length(Tp))

    # begin main solution loop
    while !terminated(d)
  
        j = mod(d.status.current_iteration, length(d.I)) + 1
        m_sp, s_sp, v_p_sp, b1_sp, b2_sp = create_solve!(Subproblem(), d, d.I[j], d.λ)
        d.status.solve_time += solve_time(m_sp)
  
        if valid_solve(Subproblem(), m_sp)
            
            save_solution!(Subproblem(), d, m_sp, s_sp, v_p_sp, b1_sp, b2_sp, d.I[j])
            update_norm_step!(d)

            if use_problem(FeasibilityProblem(), d)
                m_f = create_solve!(FeasibilityProblem(), d)
                d.status.heurestic_time += solve_time(m_f)
                if valid_solve(FeasibilityProblem(), m_f)
                    d.status.upper_bound[k] = objective_value(m_f)
                end
            end
            d.status.current_M += 1
            
            if use_problem(StepsizeProblem(), d)  # DONE
                lambda[M] .= mult    
                new_maxest = current_step(d)*current_norm(d)/alpha_step(d) + lower_bound(d)
                (d.status.maxest < new_maxest)   && (d.status.maxest = new_maxest)
                (d.status.maxest > d.status.current_estimate) && (d.status.maxest = d.status.current_estimate)
                m_ss = create_solve!(StepsizeProblem(), d, lambda, M, sstep)
                if (d.status.current_M > 5) && (d.status.current_iteration > 50) && (d.status.current_estimate > lower_bound(d))
                    if !valid_solve(StepsizeProblem(), m_ss) || (d.status.current_M >= 50000)
                        d.status.est = d.status.maxest 
                        d.status.current_step /= 10
                        d.status.current_M = 1
                        d.status.maxest = -100000
                    end
                end              
            end
            
            d.status.current_norm = d.status.prior_norm
            d.status.current_step = d.status.prior_step
            @show d.I[j]
            display_iteration(d, j)
            d.status.current_iteration += 1
        end
    end
    return nothing
end

