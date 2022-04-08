#=

=#

"""
    configure!(::Subproblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)

Set parameters used by CPLEX when solving subproblems.
"""
function configure!(::Subproblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_RANDOMSEED", j.parameter.random_seed)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", Threads.nthreads())
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 200000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    200000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    200000)
    set_time_limit_sec(m, 30)
    #if current_iteration(j) == 1
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP", 0.01)
    #end
    if current_iteration(j) > 29
        set_optimizer_attribute(m, "CPX_PARAM_INTSOLLIM", 5)
        set_optimizer_attribute(m, "CPX_PARAM_CUTUP", 0.0001)
    end
    #if current_iteration(j) == 41
     #   set_time_limit_sec(m, 30)
    #end
    return nothing
end

"""
    configure!(::FeasibilityProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)

Set parameters used by CPLEX when solving the feasibility problem.
"""
function configure!(::FeasibilityProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model)
    set_optimizer_attribute(m, "CPX_PARAM_RANDOMSEED", j.parameter.random_seed)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", Threads.nthreads())
    set_optimizer_attribute(m, "CPX_PARAM_CUTUP", current_upper_bound(j))
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP",       0.0025)
    set_optimizer_attribute(m, "CPX_PARAM_REPAIRTRIES", 1000000)
    set_optimizer_attribute(m, "CPX_PARAM_RINSHEUR",    1000000)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ",    1000000)
    return nothing
end

configure!(::StepsizeProblem, v::Val{:CPLEX}, j::JobShopProblem, m::Model) = nothing