using Distributed

# Psuedocost priority for Job i. Let `c[J]` be the cost improvement associated with solve job J and let `t[J]` be the time
# spent solving problem J then `cbar[J] = c[J]/t[J]`.

# Pseudcost group for Jobs: Let `c[i]` be the cost improvment associated 
# with solving job `i` in time J[i]. 
# The general idea is grouped pseudocost solves for parallel 

# Machine, threads, and job grouping optimization...

struct ProbRequest
    I::Vector{Int}
    flag::Bool
end

"""
Functor that produces jobs to run on workers.
"""
mutable struct SubproblemProducer
    jsp::JobShopProblem
    I::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    ni::Int
    iteration_limit::Int
    nl::Int
    nl_interval::Int
end

"""
Constructs the Producer from a JobShopProblem
"""
function SubproblemProducer(js::JobShopProblem)
end

"""
Checks termination based on producer state.
"""
terminated(p::SubproblemProducer) = terminated(p.jsp)  # DONE

"""
Creates a job (either subproblem or feasibility) until the algorithm has terminated.
"""
function (p::SubproblemProducer)(c::Channel)           # DONE
    while true
        terminated(p)   && (return nothing)
        ssp_flag = use_problem(StepsizeProblem(), p) && (p.nl != p.js.status.current_iteration)
        put!(c, ProbRequest(p.I[i], ssp_flag))
        p.ni = (p.ni > length(p.I))   ? 1 : (p.ni + 1)
    end
end

function solve_subproblem!(::Val{false}, sp::SubproblemProducer, x::ProbRequest)
    check_termination!() && return false
    m_sp, s_sp, b1_sp, b2_sp = create_solve!(Subproblem(), d, x.I, d.λ)
    d.status.solve_time += solve_time(m_sp)
    if valid_solve(Subproblem(), m_sp)
        save_solution!(Subproblem(), d, m_sp, s_sp, b1_sp, b2_sp, x.I)
        update_norm_step!(d)
        if use_problem(FeasibilityProblem(), d)
            m_f = create_solve!(FeasibilityProblem(), d)
            d.status.heurestic_time += solve_time(m_f)
            if valid_solve(FeasibilityProblem(), m_f)
                d.status.upper_bound[k] = objective_value(m_f)
            end
        end
        d.status.current_M += 1
    end
    d.status.current_norm = d.status.prior_norm
    d.status.current_step = d.status.prior_step
    display_iteration(d, j)
    return false
end
function solve_subproblem!(::Val{true}, sp::SubproblemProducer, x::ProbRequest)
    if use_problem(StepsizeProblem(), d)
        lambda[M] .= mult    
        new_maxest = current_step(d)*current_norm(d)/alpha_step(d) + lower_bound(d)
        (d.status.maxest < new_maxest)   && (d.status.maxest = new_maxest)
        (d.status.maxest > d.status.current_estimate) && (d.status.maxest = d.status.current_estimate)
        m_ss = create_solve!(StepsizeProblem(), d, lambda, M, sstep)
        if (d.status.current_M > 5) && (d.status.current_iteration > 50) && (d.status.current_estimate > lower_bound(d))
            if !valid_solve(StepsizeProblem(), m_ss) || (d.status.current_M >= 50000)
                d.status.est = d.status.maxest 
                d.status.current_step /= 10
                d.status.current_M = 1
                d.status.maxest = -100000
            end
        end              
    end
    return false
end
solve_subproblem!(sp::SubproblemProducer, x::ProbRequest) = solve_subproblem!(x.flag, sp, x)

# TODO: CHECK CONNECTION BETWEEN put! and threads
function async_solve(js::JobShopProblem)
    p = SubproblemProducer(js)
    pc = Channel(p)
    sb = Channel()
    for x in pc
        @async put!(sb, () -> solve_subproblem!(p, x))
    end
    for i in workers()
        rmprocs(i)
    end
    return
end


