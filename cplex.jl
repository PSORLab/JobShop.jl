function configure!(::Subproblem, v::Val{:CPLEX}, js, model)
    @unpack current_iteration = js.status
    set_optimizer_attribute(model, "CPX_PARAM_REPAIRTRIES", 200000)
    set_optimizer_attribute(model, "CPX_PARAM_RINSHEUR",    200000)
    set_optimizer_attribute(model, "CPX_PARAM_HEURFREQ",    200000)
    if current_iteration == 1
        set_time_limit_sec(60)
        set_optimizer_attribute(model, "CPX_PARAM_EPGAP", 0.01)
    end
    if current_iteration > 29
        set_optimizer_attribute(model, "CPX_PARAM_INTSOLLIM", 3)
        set_optimizer_attribute(model, "CPX_PARAM_CUTUP", 0.0001)
    end
    return
end

function configure!(::FeasibilityProblem, v::Val{:CPLEX}, js, model)
    @unpack current_iteration, upper_bound = js.status
    set_optimizer_attribute(model, "CPX_PARAM_CUTUP",   upper_bound)
    set_optimizer_attribute(model, "CPX_PARAM_REPAIRTRIES", 1000000)
    set_optimizer_attribute(model, "CPX_PARAM_RINSHEUR",    1000000)
    set_optimizer_attribute(model, "CPX_PARAM_HEURFREQ",    1000000)
    set_optimizer_attribute(model, "CPX_PARAM_EPGAP",        0.0025)
    return
end

function configure!(::StepsizeProblem, v::Val{:CPLEX}, js, model)
    @unpack current_iteration, upper_bound = js.status
    set_optimizer_attribute(model, "CPX_PARAM_CUTUP",   upper_bound)
    set_optimizer_attribute(model, "CPX_PARAM_REPAIRTRIES", 1000000)
    set_optimizer_attribute(model, "CPX_PARAM_RINSHEUR",    1000000)
    set_optimizer_attribute(model, "CPX_PARAM_HEURFREQ",    1000000)
    set_optimizer_attribute(model, "CPX_PARAM_EPGAP",        0.0025)
    return
end
