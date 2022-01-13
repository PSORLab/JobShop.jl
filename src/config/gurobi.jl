
const GUROBI_SEED = 12
function gurobi_config!(m::Model)
    set_optimizer_attribute(m, "Seed", GUROBI_SEED)
    set_optimizer_attribute(m, "Threads", Threads.nthreads())
    return nothing
end

function configure!(::Subproblem, v::Val{:Gurobi}, j::JobShopProblem, m::Model)
    gurobi_config!(m)
    if current_iteration(j) == 1
        set_time_limit_sec(m, 60)
        set_optimizer_attribute(m, "MIPGap", 0.01)
    end
    if current_iteration(j) > 29
        set_optimizer_attribute(m, "SolutionLimit", 3)
        set_optimizer_attribute(m, "Cutoff", 0.0001)
    end
    if current_iteration(j) == 41
        set_time_limit_sec(m, 120)
    end
    return
end

function configure!(::FeasibilityProblem, v::Val{:Gurobi}, j::JobShopProblem, m::Model)
    gurobi_config!(m)
    set_time_limit_sec(m, 60)
    set_optimizer_attribute(m, "Cutoff", current_upper_bound(j))
    set_optimizer_attribute(m, "MIPGap",       0.0025)
    return
end