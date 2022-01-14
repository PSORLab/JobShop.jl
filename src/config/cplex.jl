#=

=#

const CPLEX_SEED = 0

function cplex_config!(m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_RANDOMSEED", CPLEX_SEED)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", Threads.nthreads())
    return nothing
end

function configure!(::Subproblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    cplex_config!(m)
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 200000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    200000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    200000)
    if current_iteration(j) == 1
        set_time_limit_sec(m, 60)
        #set_optimizer_attribute(m, "CPX_PARAM_EPGAP", 0.01)
    end
    if current_iteration(j) > 29
        #set_optimizer_attribute(m, "CPX_PARAM_INTSOLLIM", 3)
        #set_optimizer_attribute(m, "CPX_PARAM_CUTUP", 0.0001)
    end
    if current_iteration(j) == 41
        set_time_limit_sec(m, 120)
    end
    return nothing
end

function configure!(::FeasibilityProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    cplex_config!(m)
    set_time_limit_sec(m, 240)
    set_optimizer_attribute(m, "CPX_PARAM_CUTUP", current_upper_bound(j))
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP",       0.0025)
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 1000000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    1000000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    1000000)
    return nothing
end