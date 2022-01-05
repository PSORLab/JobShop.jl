#DONE...
"""
$TYPEDSIGNATURES

Populates constraints and variables shared by subproblems, original problem, and feasibility problem.
"""
function create_problem(::SharedProblem, jsp::JobShopProblem, I::Vector{Int})

    @unpack J, T, Tp, R, p, ps, pr, w, d, U, O = jsp

    m = Model(jsp.parameter.optimizer)
    set_silent(m)
    @variables(m, begin
        0 <= b1[i=I, j=J[i]] <= T[end], Int                  # DONE
        0 <= c1[i=I, j=J[i]] <= T[end], Int                  # DONE
        0 <= b2[i=I, j=J[i], jᵖ=J[i], r=R] <= T[end], Int    # DONE
        0 <= c2[i=I, j=J[i], jᵖ=J[i], r=R] <= T[end], Int    # DONE
        bI1[i = I, j = J[i], t = Tp], Bin                    # DONE
        bI2[i = I, j = J[i], jᵖ = J[i], r = R, t = Tp], Bin  # DONE
    end)

    for i ∈ I, j ∈ J[i]
        if haskey(U, (i,j))
            @constraints(m, begin 
                sum((t + p[i,j])*bI1[i,j,t] for t = Tp, m = U[i,j]) - c1[i,j] == 1   # DONE                              # equation 2
                sum(bI1[i,j,t] for t = Tp, m = U[i,j]) == 1                          # DONE                             # equation 3
                sum(t*bI1[i,j,t] for t = Tp, m = U[i,j]) - b1[i,j] == 0              # DONE                              # equation 4
                [jᵖ = J[i], r = R], sum((t + p[i,j])*bI2[i,j,jᵖ,r,t] for t = Tp, m = U[i,j]) - c2[i,j,jᵖ,r] == 1  # DONE # equation 5
                [ jᵖ = J[i], r = R], sum(bI2[i,j,jᵖ,r,t] for t = Tp, m = U[i,j]) == 1                             # DONE                      # equation 6
                [jᵖ = J[i], r = R], sum(t*bI2[i,j,jᵖ,r,t] for t = Tp, m = U[i,j]) - b2[i,j,jᵖ,r] == 0   # DONE          # equation 7
            end)
        end
    end

    for i ∈ I, (k,j) ∈ enumerate(J[i][1:(end-1)])                                                               # equation 8
        jn = J[i][k+1]
        @constraints(m, begin 
            b1[i,jn] - c1[i,j] >= 1
            [jᵖ = J[i], r = R], b2[i,jn,jᵖ,r] - c2[i,j,jᵖ,r] >= 1
        end)
    end 
    
    @constraints(m, begin
        eqn_09[i = I, jᵖ = J[i]], b2[i, 1,  jᵖ, 0] - c1[i, jᵖ] >= 1  # DONE
        eqn_10[i = I, jᵖ = J[i]], b2[i, jᵖ, jᵖ, 1] - c1[i, jᵖ] >= 1  # DONE
    end)

    return m, b1, c1, b2, c2, bI1, bI2
end
create_problem(::SharedProblem, jsprob::JobShopProblem) = create_problem(SharedProblem(), jsprob, jsprob.I)

"""
$TYPEDSIGNATURES

Constructs the original (not Lagragian branch and cut) version of the Jobshop problem.
"""
function create_problem(::OriginalProblem, jsprob::JobShopProblem)
    @unpack J, Tp, R, p, ps, pr, w, d, U, O, I, Mi, γ, τ = jsp

    # add constraints shared by multiple problems 
    model, b1, c1, b2, c2, bI1, bI2 = create_problem(SharedProblem(), jsp)

    @expression(model, ex_γ[i = I, j = J[i], m = M, t = Tp], sum(b1[q,j,k]   for q = (t-p[i,j,m,1]):t))
    @expression(model, ex_τ[i = I, j = J[i], m = M, t = Tp], sum(b2[q,j,k,2] for q = (t-p[i,j,m,1]):t))
    @constraint(model, [m = Mi, t = Tp], sum(γ[i,j]*ex_γ[i,j,m,t] + τ[i,j]*ex_τ[i,j,m,t] for (i,j) in O[m]) <= M[m])
    
    # add objective
    @objective(model, Min, sum(o[i] for i = I))
    return model, b1, c1, b2, c2, bI1, bI2
end

"""
$TYPEDSIGNATURES

Constructs the feasibility problem called after convergence using estimates for 
starting times `sb1` obtain by Lagrangian relaxation.
"""
function create_solve!(::FeasibilityProblem, d::JobShopProblem)
    @unpack J, w, d, U, O, I, sb1 = d
    m, b1, c1, b2, c2, bI1, bI2 = create_problem(OriginalProblem(),d)
    @constraint(m, [i = I, j = J[i]], -1 <= b1[i,j] - sb1[i,j] <= 1)
    optimize!(m)
    finalize!(backend(model))
    GC.gc()
    return m
end

function use_problem(::FeasibilityProblem, d::JobShopProblem)
    @unpack current_norm, current_iteration = d.status
    @unpack feasible_norm_limit, feasible_interval, feasible_start, verbosity = d.parameter
    flag = current_iteration >= feasible_start
    flag &= iszero(mod(current_iteration, feasible_interval))
    flag &= current_norm < feasible_norm_limit
    (flag && verbosity > 1) && println("Stepsize problem will be solved.") : println("Stepsize problem skipped.")
    flag
end