module JobShop

using CSV, DataFrames, DocStringExtensions, JuMP, PiecewiseLinearOpt, Printf, Requires, UnPack

export configure!, load_from_csv, solve

include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "load_csv.jl"))

function configure!(t::AbstractLagrangianSubproblem, v, js, m)
    @warn "No Jobshop.configure!(t::$(typeof(t)), js, m::$(typeof(m))) set. Using default solver parameters." 
end
configure!(t::AbstractLagrangianSubproblem, js, m) = configure!(t, Val(Symbol(solver_name(m))), js, m)

include(joinpath(@__DIR__, "display.jl"))

const VALID_TSTATUS = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SOLUTION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.MEMORY_LIMIT)
const VALID_PSTATUS = (MOI.FEASIBLE_POINT)

function valid_solve(::AbstractLagrangianSubproblem, m)
    if !(termination_status(m) ∈ VALID_TSTATUS) || (primal_status(m) != MOI.FEASIBLE_POINT)
        return false
    end
    return true
end

include(joinpath(@__DIR__, "feasibility_problem.jl"))
include(joinpath(@__DIR__, "subproblem.jl"))
include(joinpath(@__DIR__, "stepsize_problem.jl"))
include(joinpath(@__DIR__, "serial_solve.jl"))
#include(joinpath(@__DIR__, "parallel_solve.jl"))

function __init__()
    @require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0"      include(joinpath(@__DIR__, "config", "cplex.jl"))
    @require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b"     include(joinpath(@__DIR__, "config", "gurobi.jl"))
    #@require KNITRO="67920dd8-b58e-52a8-8622-53c4cffbe346"     include(joinpath(@__DIR__, "config", "knitro.jl"))
    #@require MosekTools="1ec41992-ff65-5c91-ac43-2df89e9693a4" include(joinpath(@__DIR__, "config", "mosek.jl"))
    #@require Xpress="9e70acf3-d6c9-5be6-b5bd-4e2c73e3e054"     include(joinpath(@__DIR__, "config", "xpress.jl"))
end

end
