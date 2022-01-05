

"""
$TYPEDSIGNATURES

Creates subproblem used to determine step-size
"""
function create_solve!(t::StepsizeProblem, d::JobShopProblem)
    @unpack λp, step_size, Mi, Tp = d
    @unpack feasible_lambda_max, feasible_window, feasible_interval = jsp.parameter
    model = Model()
    @variable(model, 0 <= λ[m = Mi, t = Tp] <= feasible_lambda_max)
    ml = size(λp,3)
    mli(i) = ml - feasible_interval*i
    for k in [mli(i) for i = 1:feasible_window if 0 < mli(i) < ml - feasible_interval]
        kn = k + feasible_interval
        c = (1 - 4*d.status.current_step)^feasible_interval
        @constraint(model, [m = Mi, t = Tp], c*(λ[m, t] - λp[m, t, k])^2 >= (λ[m, t] - λp[m, t, kn])^2)
    end
    optimize!(model)
    finalize!(backend(model))
    GC.gc()
    return model
end
function use_problem(::StepsizeProblem, d::JobShopProblem)
    flag = iszero(mod(d.status.current_iteration, d.parameter.feasible_interval))
    flag &= d.status.current_iteration > d.parameter.feasible_start
    (flag && verbosity > 1) && println("Stepsize problem will be solved.") : println("Stepsize problem skipped.")
    return flag
end
