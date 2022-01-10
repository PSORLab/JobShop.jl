function configure!(t::AbstractLagrangianSubproblem, v, j, m)
    @warn "No Jobshop.configure!(t::$(typeof(t)), js, m::$(typeof(m))) set. Using default solver parameters." 
end
configure!(t::AbstractLagrangianSubproblem, j::JobShopProblem, m::Model) = configure!(t, Val(Symbol(solver_name(m))), j, m)

const VALID_TSTATUS = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SOLUTION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.MEMORY_LIMIT)
const VALID_PSTATUS = (MOI.FEASIBLE_POINT)

function valid_solve(::AbstractLagrangianSubproblem, m)
    if !(termination_status(m) ∈ VALID_TSTATUS) || (primal_status(m) != MOI.FEASIBLE_POINT)
        return false
    end
    return true
end

function close_problem!(m::Model)
    finalize!(backend(m))
    GC.gc()
    return
end