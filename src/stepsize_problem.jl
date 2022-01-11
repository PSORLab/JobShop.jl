"""
$TYPEDSIGNATURES

Creates subproblem used to determine step-size.
"""
function solve_problem(::StepsizeProblem, j::JobShopProblem)

    stepsize_start = time()

    @unpack lambd, MachineType, T = j
    @unpack current_step, current_norm, current_M = j.status
    @unpack stepsize_lambda_max, stepsize_interval, alpha_step_2 = j.parameter
    
    new_maxest = current_step*current_norm/alpha_step_2 + current_lower_bound(j)
    (j.status.maxest < new_maxest)        && (j.status.maxest = new_maxest)
    (j.status.maxest > j.status.estimate) && (j.status.maxest = j.status.estimate)
    
    model = direct_model(optimizer_with_attributes(j.parameter.optimizer))
    configure!(StepsizeProblem(), j, model)
    set_silent(model)

    @variable(model, 0 <= λ[MachineType,T] <= stepsize_lambda_max)
    todo_param = 2
    c = (1 - 2*todo_param*current_step)^stepsize_interval
    for n = 1:100
        for k = (current_M-2-1000):(current_M-stepsize_interval)
            if k == stepsize_interval*n
                kn = k + stepsize_interval
                @constraint(model, [m=MachineType,t=T], c*(λ[m,t] - lambd[k][m,t])^2 >= (λ[m,t] - lambd[kn][m,t])^2)
            end
        end
    end
    optimize!(model)
    j.status.time_solve_stepsize += solve_time(model)
    valid_flag = valid_solve(StepsizeProblem(), model) 

    close_problem!(model)
    j.status.time_total_stepsize = time() - stepsize_start
    return valid_flag
end

function use_problem(::StepsizeProblem, j::JobShopProblem)
    @unpack stepsize_interval, stepsize_start, verbosity = j.parameter
    flag = iszero(mod(current_iteration(j), stepsize_interval))
    flag &= current_iteration(j) > stepsize_start
    (verbosity > 1) && (flag && println("Stepsize problem will be solved."))
    return flag
end