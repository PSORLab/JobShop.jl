#DONE...
"""
$TYPEDSIGNATURES

Populates constraints and variables shared by subproblems, original problem, and feasibility problem.
"""
function create_problem(::SharedProblem, jsp::JobShopProblem, I::Vector{Int})

    @unpack J, T, R, p, ps, pr, w, d, U, O = jsp

    m = Model()

    @variable(m, 0 <= b1[i ∈ I, j ∈ J[i]] <= T[end], Int)
    @variable(m, 0 <= c1[i ∈ I, j ∈ J[i]] <= T[end], Int)

    @variable(m, 0 <= b2[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], r ∈ R] <= T[end], Int)
    @variable(m, 0 <= c2[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], r ∈ R] <= T[end], Int)

    @variable(m, bI1[i ∈ I, j ∈ J[i], t ∈ T], Bin)
    @variable(m, bI2[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], t ∈ T, r ∈ R], Bin)

    @constraints(m, begin 
        eqn2[i ∈ I, j ∈ J[i]], sum((t + p[i,j])*bI1[i,j,t] for t ∈ T, m ∈ U[i,j]) - c1[i,j] == 1 
        eqn3[i ∈ I, j ∈ J[i]], sum(bI1[i,j,t] for t ∈ T, m ∈ U[i,j]) == 1
        eqn4[i ∈ I, j ∈ J[i]], sum(t*bI1[i,j,t] for t ∈ T, m ∈ U[i,j]) - b1[i,j] == 0
        eqn5[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], r ∈ R], sum((t + p[i,j])*bI2[i,j,jᵖ,t,r] for t ∈ T, m ∈ U[i,j]) - c2[i,j,jᵖ,r] == 1 
        eqn6[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], r ∈ R], sum(bI2[i,j,jᵖ,t,r] for t ∈ T, m ∈ U[i,j]) == 1
        eqn7[i ∈ I, j ∈ J[i], jᵖ ∈ J[i], r ∈ R], sum(t*bI2[i,j,jᵖ,t,r] for t ∈ T, m ∈ U[i,j]) - b2[i,j,jᵖ,r] == 0
    end)

    # equation 8
    for i ∈ I
        for (k,j) in enumerate(J[i][1:(end-1)])
            jn = J[i][k+1]
            @constraints(m, begin 
                b1[i,jn] - c1[i,j] >= 1
                [jᵖ ∈ J[i], r ∈ R], b2[i,jn,jᵖ,r] - c2[i,j,jᵖ,r] >= 1
            end)
        end
    end 
    
    @constraints(m, begin
        eqn_09[i ∈ I, jᵖ ∈ J[i]], b2[i, 1,  jᵖ, 0] - c1[i, jᵖ] >= 1
        eqn_10[i ∈ I, jᵖ ∈ J[i]], b2[i, jᵖ, jᵖ, 1] - c1[i, jᵖ] >= 1
    end)

    return m, b1, c1, b2, c2, bI1, bI2
end
create_problem(::SharedProblem, jsprob::JobShopProblem) = create_problem(SharedProblem(), jsprob, jsprob.I)

"""
$TYPEDSIGNATURES

Creates the expression used by the objective (via SOS constraint). 
"""
function original_objective_ex!(m::Model, jsprob::JobShopProblem, c1, c2, I)
    @unpack J, d, ps, pr, w, T = jsprob

    # objective
    ps_eq11a = Dict{Tuple{Int,Int},Float64}()
    t = Dict{Tuple{Int,Int},Float64}()
    ta = Dict{Int,Float64}()
    tb = Dict{Int,Float64}()
    tc = Dict{Int,Float64}()
    for i ∈ I
        jm1 = J[i][1:(end-1)]
        for j ∈ J[i]
            t[i,j] = prod(k -> (1 - ps[i,k]), jm1) - prod(k -> (1 - ps[i,k]), J[i])
        end
        ta[i] = prod(j -> (1 - ps[i,j]),             J[i])
        tb[i] = sum(j -> t[i,j]*(1 - pr[i,j]), J[i])
        tc[i] = sum(j -> t[i,j]*pr[i,j],       J[i])
    end
    @variable(m, α[i ∈ I,1:3])
    @constraints(m, begin 
        [i ∈ I], α[i,1] - c1[i,J[i][end]]   + d[i] == 0                # TODO: check how to match up c1[i,j,j',r] index to completed value...
        [i ∈ I], α[i,2] - c2[i,J[i][end], J[i][end-1], 0] + d[i] == 0
        [i ∈ I], α[i,3] - c2[i,J[i][end], J[i][end-1], 1] + d[i] == 0
    end)
    za = Dict{Int,Any}()
    zb = Dict{Int,Any}()
    zc = Dict{Int,Any}()
    for i ∈ I
        za[i] = piecewiselinear(m, α[i,1], [-T[end],0,2*T[end]], u -> max(u, 0.0))
        zb[i] = piecewiselinear(m, α[i,2], [-T[end],0,2*T[end]], u -> max(u, 0.0))
        zc[i] = piecewiselinear(m, α[i,3], [-T[end],0,2*T[end]], u -> max(u, 0.0))
    end
    @expression(m, o[i ∈ I], w[i]*ta[i]*za[i] + w[i]*tb[i]*zb[i] + w[i]*tc[i]*zc[i])
    return o
end

"""
$TYPEDSIGNATURES

Constructs the original (not Lagragian branch and cut) version of the Jobshop problem.
"""
function create_problem(::OriginalProblem, jsprob::JobShopProblem)
    @unpack J, T, R, p, ps, pr, w, d, U, O, I, Mi = jsp

    # add constraints shared by multiple problems 
    model, b1, c1, b2, c2, bI1, bI2 = create_problem(SharedProblem(), jsp)

    # add machine capacity constraints
    γ = Dict{Tuple{Int,Int},Float64}()
    τ = Dict{Tuple{Int,Int},Float64}()
    for i ∈ I, j ∈ J[i]
        jm1 = J[i][1:(end-1)]
        if (i,j) in O[m]
            γ[i,j] = prod(k -> 1 - ps[i,k], jm1)
            τ[i,j] = prod(k -> 1 - ps[i,k], jm1) - prod(k -> 1 - ps[i,k], J[i])*pr[i,j]
        end
    end
    @expression(model, ex_γ[i ∈ I, j ∈ J[i], m ∈ M, t ∈ T], sum(b1[q,j,k]   for q ∈ (t-p[i,j,m,1]):t))
    @expression(model, ex_τ[i ∈ I, j ∈ J[i], m ∈ M, t ∈ T], sum(b2[q,j,k,2] for q ∈ (t-p[i,j,m,1]):t))
    @constraint(model, [m ∈ Mi, t ∈ T], sum(γ[i,j]*ex_γ[i,j,m,t] + τ[i,j]*ex_τ[i,j,m,t] for (i,j) in O[m]) <= M[m])
    
    # add objective
    o = original_objective_ex!(model, jsp, c1, c2, jsp.I)
    @objective(model, Min, sum(o[i] for i ∈ I))
    return model, b1, c1, b2, c2, bI1, bI2
end

# DONE.
"""
$TYPEDSIGNATURES

Constructs the feasibility problem called after convergence using estimates for 
starting times `sb1` obtain by Lagrangian relaxation.
"""
function create_problem(::FeasibilityProblem, d::JobShopProblem)
    @unpack J, T, R, p, ps, pr, w, d, U, O, I = d
    model, b1, c1, b2, c2, bI1, bI2 = create_problem(OriginalProblem(),d)
    d.bpm1 = @constraint(model, bpm1[i ∈ I, j ∈ J[i]], -1 <= b1[i,j] - 0.0 <= 1)
    d.feasibility_model = model
    return nothing
end

function update_problem(::FeasibilityProblem, d::JobShopProblem)
    @unpack feasibility_model, sb1, bpm1 = d
    delete(feasibility_model, bpm1)
    d.bpm1 = @constraint(feasibility_model, bpm1[i ∈ I, j ∈ J[i]], -1 <= b1[i,j] - sb1[i,j] <= 1)
    return nothing
end

function update_solve!(::FeasibilityProblem, d::JobShopProblem)
    update_subproblem!(FeasibilityProblem(), d)
    optimize!(d.feasibility_model)
    return solve_time(d.feasibility_model)
end

function use_problem(::FeasibilityProblem, d::JobShopProblem)
    check_feasible_lambda(d.status.currrent_iteration, d.parameter.feasible_interval)
end