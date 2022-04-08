

abstract type AbstractExt end

struct Ext <: AbstractExt end

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

Type used to indicate the stepsize problem solved.
"""
struct StepsizeProblem <: AbstractLagrangianSubproblem end

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

"""
$(TYPEDEF)

Parameters used in Lagrangian solution method.
"""
Base.@kwdef mutable struct SolveParameter
    "Seed used by Random package. Set as a parameter for reproducibility."
    random_seed::Int = 0
    "The number of time steps to search about the Lagragian dual solution for feasible points."
    feasibility_window::Int = 3
    "Starting norm."
    start_norm::Float64 = 100.0
    "Starting stepsize."
    start_step::Float64 = 200.0
    "Probability of scrap."
    prob::Float64           = 0.05
    "Probability of rework."
    prob_r::Float64         = 0.2
    "Shift length."
    ShiftLength::Int        = 18
    "Starting penalty."
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
    "Absolute threshold for the norm below which a feasible problem solution may be attempted."
    feasible_norm_limit::Float64 = 10.0
    "Number of iterations wait before checking to see if a feasible problem should be solved."
    feasible_start::Int = 100
    "Number of times to restart the feasible problem solution"
    feasible_solve_count::Int = 3
    "Amount of time spent to solve a feasible problem."
    feasible_solve_time::Float64 = 240.0
    "Printing verbosity."
    verbosity::Int = 1
    "Optimizer to used in Lagrangian relaxations."
    optimizer = nothing
    "Number of iterations between updating the penalty term."
    penalty_iteration::Int = 200
    "Factor used to multiply the penalty by when updating the penalty term."
    penalty_factor::Float64 = 1.05
    "Stepsize problem set by solving optimization problem"
    use_stepsize_program::Bool = false
    "Subproblems to solve prior to computing new stepsize"
    stepsize_interval::Int = 20
    alpha_step::Float64 = 0.5
    ext::AbstractExt = Ext()
end

"""
$(TYPEDEF)

Storage to characterize the state of the solution routine.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolveStatus
    "Current penalty."
    penalty::Float64 = 120.0
    "Iteration number."
    current_iteration::Int = 0
    "Current norm"
    current_norm::Float64 = 100.0
    "Current step size."
    current_step::Float64 = 0.1
    "Norm of previous iteration."
    prior_norm::Float64 = 100.0
    "Step size from previous iteration."
    prior_step::Float64 = 0.1
    "Stores lower bound at each iteration"
    lower_bound::Dict{Int,Float64} = Dict{Int,Float64}()
    "Stores upper bound at each iteration"
    upper_bound::Dict{Int,Float64} = Dict{Int,Float64}()
    "Total lower bounding calculation time by iteration."
    lower_bound_time::Dict{Int,Float64} = Dict{Int,Float64}()
    "Total upper bounding calculation time by iteration."
    upper_bound_time::Dict{Int,Float64} = Dict{Int,Float64}()
    "Current estimate."
    estimate::Float64 = 0.0
    "Current maximum of the estimate"
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
    "Flag to indicate that a feasible solution was found."
    feasible_problem_found::Bool = false
    "Contains current M value used in stepsize calculation."
    current_M::Int = 1
end

Base.@kwdef mutable struct JobShopProblem
    "Maximum time window considered."
    Tmax::Int = 1000
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
    lambd::Dict{Int,Matrix{Float64}} = Dict{Int,Matrix{Float64}}()
    sTard1::Vector{Float64} = Float64[]
    sTard2::Array{Float64,3} = zeros(Float64,2,2,2)
    T::UnitRange{Int} = 1:232
    sbTime1::Matrix{Float64}  = zeros(Float64,2,2)
    sbTime2::Dict{Tuple{Int,Int,Int,Int},Float64} = Dict{Tuple{Int,Int,Int,Int},Float64}()
    status::SolveStatus       = SolveStatus()
    parameter::SolveParameter = SolveParameter()
    sbI1_indx::Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}} = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}}()
    sbI2_indx::Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}} = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}}()
    sbI1_indx_n::Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}} = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}}()
    sbI2_indx_n::Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}} = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}}()
    feasible_bTime1::Dict{Any,Int} = Dict{Any,Int}()
    feasible_bTime2::Dict{Any,Int} = Dict{Any,Int}()
end

function initialize!(j::JobShopProblem)
    @unpack MachineType, T, I, Ii, IJT, MIJ, Jop = j

    j.mult = zeros(length(MachineType), length(T))
    Jmax = maximum(x -> x.j, IJT)
    j.sbTime1 = zeros(Float64,length(I), Jmax)
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

    sbI1_indx = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}}()
    sbI2_indx = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}}()
    sbI1_indx_n = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int}}}}()
    sbI2_indx_n = Dict{Int,Matrix{Vector{Tuple{Int,Int,Int,Int}}}}()
    for i = 1:length(Ii)
        sbI1_indx[i] = Array{Vector{Tuple{Int,Int,Int}}}(undef, length(MachineType), length(T))
        sbI2_indx[i] = Array{Vector{Tuple{Int,Int,Int,Int}}}(undef, length(MachineType), length(T))
        sbI1_indx_n[i] = Array{Vector{Tuple{Int,Int,Int}}}(undef, length(MachineType), length(T))
        sbI2_indx_n[i] = Array{Vector{Tuple{Int,Int,Int,Int}}}(undef, length(MachineType), length(T))
        for m in MachineType, t in T
            sbI1_indx[i][m,t] = Tuple{Int,Int,Int}[]
            sbI2_indx[i][m,t] = Tuple{Int,Int,Int,Int}[]
            sbI1_indx_n[i][m,t] = Tuple{Int,Int,Int}[]
            sbI2_indx_n[i][m,t] = Tuple{Int,Int,Int,Int}[]
        end
        outi = fill(true, length(I))
        for k in Ii[i]
            outi[k] = false
        end
        for mi in MachineType, t in T, Ma in MIJ, p in IJT
            if (p.i == Ma.i) && (p.j == Ma.j) && (Ma.m == mi)
                if outi[p.i]
                    for k in (t - p.t + 1):t
                        if t - p.t >= 0
                            push!(sbI1_indx[i][mi,t], (Ma.i,Ma.j,k))
                        end
                    end
                else
                    for k in (t - p.t + 1):t
                        if t - p.t >= 0
                            push!(sbI1_indx_n[i][mi,t], (Ma.i,Ma.j,k))
                        end
                    end
                end
            end
            for j in Jop[Ma.i]
                if (p.i == Ma.i) && (p.j == Ma.j) && (Ma.m == mi)
                    if outi[p.i]
                        for k in (t - p.t+1):t
                            if t - p.t >= 0
                                push!(sbI2_indx[i][mi,t], (Ma.i,Ma.j,j,k))
                            end
                        end
                    else
                        for k in (t - p.t+1):t
                            if t - p.t >= 0
                                push!(sbI2_indx_n[i][mi,t], (Ma.i,Ma.j,j,k))
                            end
                        end
                    end
                end
            end
        end
    end
    j.sbI1_indx = sbI1_indx
    j.sbI2_indx = sbI2_indx
    j.sbI1_indx_n = sbI1_indx_n
    j.sbI2_indx_n = sbI2_indx_n
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