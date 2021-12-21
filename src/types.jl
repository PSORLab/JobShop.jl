


"""
$(TYPEDEF)

Abstract type used to describe each subproblem
"""
abstract type AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate components shared by original, feasibility, and subproblems.
"""
struct SharedProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the original problem.
"""
struct OriginalProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the feasibility problem solved at termination.
"""
struct FeasibilityProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the problem solved to determine stepsize.
"""
struct StepsizeProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the subproblem.
"""
struct Subproblem <: AbstractLagrangianSubproblem end

"""
"""
struct CoordinatingProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

User specified parameters set when attempting to solve the problem.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolveParameter
    "Total iteration limit"
    iteration_limit::Int = 60000
    "If using a sequential solution routine, store subproblem  model and incrementally update if `true`"
    store_subproblems::Bool = false
    "Upper bound for dual values used in feasibility problem formulation"
    feasible_lambda_max = 232.0  #TODO: Why is this the maximal value?
    "Absolute tolerance criteria for termination"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance criteria for termination"
    relative_tolerance::Float64 = 1E-3
end

"""
$(TYPEDEF)

Storage to characterize the state of the solution routine.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolveStatus
    "Total time spent solving the problem"
    total_solve_time = 0
    "Current iteration the problem has been solved to"
    current_iteration::Int = 0
    "Stores lower bound at each iteration"
    lower_bound::Dict{Int,Float64} = Dict{Int,Float64}()
    "Stores upper bound at each iteration"
    upper_bound::Dict{Int,Float64} = Dict{Int,Float64}()
end

"""
$(TYPEDEF)

Main structure used to hold problem specifications and solutions status.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct JobShopProblem
    "Process time"
    p::Dict{Tuple{Int,Int},Int} = Dict{Tuple{Int,Int},Int}()
    "List of operations involving for each I"
    J::Dict{Int,Vector{Int}}        = Dict{Int,Vector{Int}}()
    "Time Range"
    T::UnitRange{Int}               = 1:1000
    "Rework statuses"
    R::UnitRange{Int}               = 0:1
    "Machine capacity"
    M::Vector{Int}                  = Int[]
    "Time dart is due"
    d::Vector{Int}                  = Int[]
    "Machines that can run (part,operation) combinations"
    U::Dict{Tuple{Int,Int},Vector{Int}} = Dict{Tuple{Int,Int},Vector{Int}}()
    "Allowable (part operation) combinations per machine"
    O::Dict{Int,Vector{Tuple{Int,Int}}} = Dict{Int,Vector{Tuple{Int,Int}}}()
    o::Dict{Int, Any}                   = Dict{Int, Any}()
    g::Dict{Int, Any}                   = Dict{Int, Any}()
    s::Dict{Int, Any}                   = Dict{Int, Any}()
    "Subproblem Storage (if used)"
    m::Dict{Int, Any}               = Dict{Int, Any}()
    "Parameters used in solving jobshop problem"
    parameter::SolveParameter       = SolveParameter()
    "Status of jobshop problem solution"
    status::SolveStatus             = SolveStatus()
end

lower_bound(d::JobShopProblem) = d.status.lower_bound[d.status.current_iteration]
upper_bound(d::JobShopProblem) = d.status.upper_bound[d.status.current_iteration]
current_abs_tol(d::JobShopProblem) = upper_bound(d) - lower_bound(d)
function current_rel_tol(d::JobShopProblem)
    L = lower_bound(d)
    U = upper_bound(d)
    return abs(U - L)/(max(abs(L), abs(U))) : Inf
end
