"""
$TYPEDSIGNATURES

Creates subproblem used to determine step-size.
"""
function solve_problem(::StepsizeProblem, jsp::JobShopProblem)

    stepsize_start = time()

    @unpack MachineType, T, Tmax, lambd = jsp
    @unpack current_step, current_norm = jsp.status
    @unpack stepsize_interval = jsp.parameter
    
    jsp.status.current_M += 1
    current_M = jsp.status.current_M

    model = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    configure!(StepsizeProblem(), jsp, model)
    set_silent(model)

    # save multiplier
    for m=MachineType,t=T
        lambd[current_M] = copy(jsp.mult)
    end 

    @variable(model, 0 <= λ[MachineType,T] <= Tmax)
    todo_param = 2
    c = (1 - 2*todo_param*current_step)^stepsize_interval
    for n = 1:100
        for k = (current_M-2-1000):(current_M-stepsize_interval)
            if k == stepsize_interval*n
                kn = k + stepsize_interval
                @constraint(model, c*sum((λ[m,t] - lambd[k][m,t])^2 for m=MachineType, t=T) >= sum((λ[m,t] - lambd[kn][m,t])^2 for m=MachineType,t=T))
            end
        end
    end
    optimize!(model)
    jsp.status.time_solve_stepsize += solve_time(model)
    valid_flag = valid_solve(StepsizeProblem(), model) 

    if jsp.status.current_iteration > 25 
        if jsp.status.maxest < jsp.status.current_step*jsp.status.current_norm/jsp.parameter.alpha_step + current_lower_bound(jsp)
            jsp.status.maxest = jsp.status.current_step*jsp.status.current_norm/jsp.parameter.alpha_step + current_lower_bound(jsp)
        end
    end
    if jsp.status.maxest > jsp.status.estimate
        jsp.status.maxest = jsp.status.estimate
    end
    if (jsp.status.current_M > 5) && (jsp.status.current_iteration > 50) && (jsp.status.estimate > current_lower_bound(jsp))
        if !valid_flag || (jsp.status.current_M >= 50000)
            jsp.status.estimate = jsp.status.maxest 
            jsp.status.current_M = 1
            jsp.status.maxest = -100000
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
    (verbosity > 1) && (flag && println("Stepsize problem will be solved."))
    return flag
end