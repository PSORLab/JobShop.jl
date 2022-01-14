
"""
$(TYPEDEF)

Abstract type used to describe each subproblem
"""
abstract type AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the feasibility problem solved at termination.
"""
struct FeasibilityProblem <: AbstractLagrangianSubproblem end

"""
$(TYPEDEF)

Type used to indicate the subproblem.
"""
struct Subproblem <: AbstractLagrangianSubproblem end

struct PartOpT
    i::Int
    j::Int
    t::Int
end
has_ij(x::PartOpT, i, j) = x.i == i && x.j == j

struct MaPartOp
    m::Int
    i::Int
    j::Int
end
has_nmi(P::PartOpT, Ma::MaPartOp, mi, Ii) = P.i == Ma.i && P.j == Ma.j && Ma.m == mi && !(P.i in Ii)
has_mi(P::PartOpT, Ma::MaPartOp, mi, Ii)  = P.i == Ma.i && P.j == Ma.j && Ma.m == mi && (P.i in Ii)

Base.@kwdef mutable struct SolveParameter
    feasibility_window::Int = 3
    start_norm::Float64 = 100.0
    start_step::Float64 = 200.0
    prob::Float64           = 0.05
    prob_r::Float64         = 0.2
    ShiftLength::Int        = 18
    penalty::Float64 = 120.0
    "starting estimate"
    estimate::Float64 = 1400.0
    "Absolute tolerance criteria for termination"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance criteria for termination"
    relative_tolerance::Float64 = 1E-3
    "Total iteration limit"
    iteration_limit::Int = 5000
    "Starting upper bound"
    start_upper_bound::Float64 = Inf
    feasible_norm_limit::Float64 = 10.0
    "feasible lambda"
    feasible_start::Int = 100
    verbosity::Int = 1
    optimizer = nothing
    penalty_increase_iteration::Int = 4000
end

"""
$(TYPEDEF)

Storage to characterize the state of the solution routine.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolveStatus
    penalty::Float64 = 120.0
    current_iteration::Int = 0
    current_norm::Float64 = 100.0
    current_step::Float64 = 0.006
    prior_norm::Float64 = 100.0
    prior_step::Float64 = 0.006
    "Stores lower bound at each iteration"
    lower_bound::Dict{Int,Float64} = Dict{Int,Float64}()
    "Stores upper bound at each iteration"
    upper_bound::Dict{Int,Float64} = Dict{Int,Float64}()
    lower_bound_time::Dict{Int,Float64} = Dict{Int,Float64}()
    upper_bound_time::Dict{Int,Float64} = Dict{Int,Float64}()
    estimate::Float64 = 0.0
    maxest::Float64 = -100000
    "Time at which solution algorithm begins"
    time_start::Float64 = 0.0
    "Time spent by optimizer solving subproblems"
    time_solve_subprob::Float64 = 0.0
    "Total time spent solving subproblems"
    time_total_subprob::Float64 = 0.0
    "Time spent by optimizer solving stepsize problems"
    time_solve_stepsize::Float64 = 0.0
    "Total time spent solving stepsize problems"
    time_total_stepsize::Float64 = 0.0
    "Time spent by optimizer solving feasibility problems"
    time_solve_feasibility::Float64 = 0.0
    "Total time spent solving feasibility problems"
    time_total_feasibility::Float64 = 0.0
    feasible_problem_found::Bool = false
end

Base.@kwdef mutable struct JobShopProblem
    "Time dart is due"
    PartDue::Vector{Int}    = Int[]
    MachineCap::Vector{Int} = Int[]
    MachineType::UnitRange{Int} = 1:1
    R::UnitRange{Int} = 1:2
    I::Vector{Int} = Int[]
    Ii::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    J::Dict{Int,Int} = Dict{Int,Int}()
    Jop::Dict{Int,UnitRange{Int}} = Dict{Int,UnitRange{Int}}()
    IJT::Vector{PartOpT}  = PartOpT[]
    MIJ::Vector{MaPartOp} = MaPartOp[]
    "sbTimeI1"
    sbI1::Dict{Tuple{Int,Int,Int},Float64} = Dict{Tuple{Int,Int,Int},Float64}()
    "sbTimeI2"
    sbI2::Dict{Tuple{Int,Int,Int,Int,Int},Float64} = Dict{Tuple{Int,Int,Int,Int,Int},Float64}()
    sslackk::Matrix{Float64} = zeros(Float64,2,2)
    sv_p::Matrix{Float64} = zeros(Float64,2,2)
    mult::Matrix{Float64} = zeros(Float64,2,2)
    sTard1::Vector{Float64} = Float64[]
    sTard2::Array{Float64,3} = zeros(Float64,2,2,2)
    T::UnitRange{Int} = 1:220
    sbTime1 = nothing
    status::SolveStatus       = SolveStatus()
    parameter::SolveParameter = SolveParameter()
end

function initialize!(j::JobShopProblem)
    @unpack MachineType, T = j
    j.mult = zeros(length(MachineType), length(T))
    j.status.current_iteration = 1
    j.status.lower_bound[0] = -Inf
    j.status.upper_bound[0] = j.parameter.start_upper_bound
    j.status.estimate = j.parameter.estimate
    j.status.current_norm = j.parameter.start_norm
    j.status.current_step = j.parameter.start_step
    j.status.prior_norm = j.parameter.start_norm
    j.status.prior_step = j.parameter.start_step
    j.status.penalty = j.parameter.penalty
    j.status.time_start = time()
    return
end

current_iteration(d::JobShopProblem) = d.status.current_iteration
function current_lower_bound(d::JobShopProblem)
    @unpack lower_bound = d.status
    return lower_bound[maximum(keys(lower_bound))]
end
function current_upper_bound(d::JobShopProblem)
    @unpack upper_bound = d.status
    return upper_bound[maximum(keys(upper_bound))]
end
function current_abs_gap(d::JobShopProblem)
    gap = current_upper_bound(d) - current_lower_bound(d)
    return isnan(gap) ? Inf : gap
end
function current_rel_gap(d::JobShopProblem)
    L = current_lower_bound(d)
    U = current_upper_bound(d)
    if isnan(U - L) || iszero(max(abs(L), abs(U)))
        return Inf
    end
    return abs(U - L)/(max(abs(L), abs(U)))
end