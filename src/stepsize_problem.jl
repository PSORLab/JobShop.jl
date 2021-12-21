

"""
$TYPEDSIGNATURES

Creates subproblem used to determine step-size
"""
function create_problem(t::StepsizeProblem, d::JobShopProblem)
    @unpack λp, step_size =  d
    @unpack feasible_lambda_max, feasible_window, feasible_interval = jsp.parameter
    model = Model()
    @variable(model, 0 <= λ[m ∈ M, t ∈ T] <= feasible_lambda_max)
    ml = size(λp,3)
    mli(i) = ml - feasible_interval*i
    for k in [mli(i) for i = 1:feasible_window if 0 < mli(i) < ml - feasible_interval]
        kn = k + feasible_interval
        c = (1 - 4*d.status.current_step)^feasible_interval
        @constraint(model, [m ∈ M, t ∈ T], c*(λ[m, T] - λp[m, T, k])^2 >= (λ[m, T] - λp[m, T, kn])^2)
    end
    return model
end
function use_problem(::StepsizeProblem, d::JobShopProblem)
    flag = iszero(mod(d.status.currrent_iteration, d.parameter.feasible_interval))
    flag &= d.status.currrent_iteration > d.parameter.start_feasible
    return flag
end
