"""
$TYPEDSIGNATURES

Creates subproblem used to determine step-size.
"""
function solve_problem(::StepsizeProblem, j::JobShopProblem)

    stepsize_start = time()

    @unpack lambd, MachineType, T = j
    @unpack current_step, current_norm, current_M = j.status
    @unpack feasible_lambda_max, stepsize_interval, alpha_step_2 = j.parameter
    
    new_maxest = current_step*current_norm/alpha_step_2 + current_lower_bound(j)
    (j.status.maxest < new_maxest)        && (j.status.maxest = new_maxest)
    (j.status.maxest > j.status.estimate) && (j.status.maxest = d.status.estimate)
    
    model = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    @variable(model, 0 <= λ[MachineType,T] <= stepsize_lambda_max)
    c = (1 - 4*current_step)^stepsize_interval
    for n = 1:100
        for m = (current_M-2-1000):(current_M-stepsize_interval)
            if m == stepsize_interval*n
                kn = k + stepsize_interval
                @constraint(model, [m=MachineType,t=Tp], c*(λ[m,t] - lambd[k][m,t])^2 >= (λ[m,t] - lambd[kn][m,t])^2)
            end
        end
    end
    optimize!(model)
    jsp.status.time_solve_stepsize += solve_time(m)
    valid_flag = valid_solve(StepsizeProblem(), model) 
    if (current_M > 5) && (current_iteration(j) > 50) && (j.status.estimate > current_lower_bound(j))
        if !valid_solve(model) || (j.status.current_M >= 50000)
            j.status.estimate = j.status.maxest 
            j.status.current_step /= 10
            j.status.current_M = 1
            j.status.maxest = -100000
        end
    end  

    close_problem!(model)
    jsp.status.time_total_stepsize = time() - stepsize_start
    return valid_flag
end

function use_problem(::StepsizeProblem, j::JobShopProblem)
    @unpack stepsize_interval, stepsize_start, verbosity = j.parameter
    flag = iszero(mod(current_iteration(j), stepsize_interval))
    flag &= current_iteration(j) > stepsize_start
    (verbosity > 1) && (flag ? println("Stepsize problem will be solved.") : println("Stepsize problem skipped."))
    return flag
end