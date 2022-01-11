#=

=#

function configure!(::Subproblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 200000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    200000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    200000)
    if current_iteration(j) == 1
        set_time_limit_sec(m, 60)
        set_optimizer_attribute(m, "CPX_PARAM_EPGAP", 0.01)
    end
    if current_iteration(j) > 29
        set_optimizer_attribute(m, "CPX_PARAM_INTSOLLIM", 3)
        set_optimizer_attribute(m, "CPX_PARAM_CUTUP", 0.0001)
    end
    if current_iteration(j) == 41
        set_time_limit_sec(m, 120)
    end
    return
end

function configure!(::FeasibilityProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_CUTUP", current_upper_bound(j))
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 1000000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    1000000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    1000000)
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP",        0.0025)
    return
end

function configure!(::StepsizeProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_CUTUP", 0.0)
    return
end
