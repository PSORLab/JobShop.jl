"""
$TYPEDSIGNATURES

A subproblem with `i ∈ I`.
"""
function create_problem(::Subproblem, jsprob::JobShopProblem, I::Vector{Int}, λ)
    @unpack J, Tp, T, R, p, ps, pr, w, d, U, Om, Mi, M, y, τ, sslackk, stard1, stard2, ta, tb, tc, sv_p = jsprob
    @unpack penalty = jsprob.status
    model, b1, c1, b2, c2, bI1, bI2 = create_problem(SharedProblem(), jsprob, I)

    @variable(model,   0 <=  v_p[m=Mi,t=Tp] <= 16)
    @variable(model, -16 <= slackk[m=Mi,t=Tp]  <= 16) # slackk
    @variable(model,   0 <= slackk1[m=Mi,t=Tp] <= 16) # slackk1

    @expression(model, ex_y[i = I, j = J[i], m = Mi, t = Tp], sum(bI1[i,j,k]     for k = (t-p[i,j]):t if t-p[i,j] ∈ Tp))
    @expression(model, ex_τ[i = I, j = J[i], m = Mi, t = Tp], sum(bI2[i,j,j,1,k] for k = (t-p[i,j]):t if t-p[i,j] ∈ Tp))
    @expression(model, g[i = I, m = Mi, t = Tp], sum(y[i,j]*ex_y[i,j,m,t] + τ[i,j]*ex_τ[i,j,m,t] for (k,j) ∈ Om if (k == i) && (t-p[i,j] ∈ Tp) && !iszero(y[i,j])))

    @constraint(model, slackk - v_p .<= 0)
    @constraint(model, slackk + v_p .>= 0)
    @constraint(model, gp_mt[m = Mi,t = Tp], sum(g[i,m,t] for i=I) - slackk[m,t] + slackk1[m,t] == M[m])


    @variable(model, 0 <= ComTime1[i = I] <= T[end])
    @variable(model, 0 <= ComTime2[i = I, j = J[i], r = R] <= T[end])

    for i = I
        jn = J[i][end]
        @constraint(model, c1[i,jn] - ComTime1[i] == 0)
        @constraint(model, [j = J[i], r = R], ComTime2[i,j,r] - c2[i,jn,j,r] == 0)
    end

    @variable(model, 0 <= W1[1:3, i = I] <= 1)
    @variable(model, α1[1:3, i = I], Bin)

    @constraint(model, W1 .- α1 .<= 0)
    @constraint(model, [i = I], W1[1,i] + W1[2,i] + W1[3,i] == 1)
    @constraint(model, [i = I], α1[1,i] + α1[3,i]        <= 1)
    @constraint(model, [i = I], α1[2,i] == 1)
    @constraint(model, [i = I], ComTime1[i] - d[i] == W1[1,i]*(-Tp[end]) + W1[3,i]*2*Tp[end])

    @variable(model, 0 <= W2[1:3, i = I, j = J[i], r = R] <= 1)
    @variable(model, α2[1:3, i = I, j = J[i], r = R])

    @constraint(model, [k = 1:3, i = I, j = J[i], r = R], W2[k,i,j,r] - α2[k,i,j,r]  <= 0)
    @constraint(model, [i = I, j = J[i], r = R], W2[1,i,j,r] + W2[2,i,j,r] + W2[3,i,j,r] == 1)
    @constraint(model, [i = I, j = J[i], r = R], α2[1,i,j,r] + α2[3,i,j,r] <= 1) 
    @constraint(model, [i = I, j = J[i], r = R], α2[2,i,j,r] == 1)
    @constraint(model, [i = I, j = J[i], r = R], ComTime2[i,j,r] - d[i] == W2[1,i,j,r]*(-Tp[end]) + W2[3,i,j,r]*2*Tp[end])

    @expression(model, tard1[i = I], W1[3,i]*2*Tp[end])
    @expression(model, tard2[i = I, j = J[i], r = R], W2[3,i,j,r]*2*Tp[end])
    @expression(model, total_tard, sum(ta[i]*tard1[i] for i = I) +
                                   sum(tb[i]*tard2[i,j,1] + tc[i]*tard2[i,j,2] for i = I, j = J[i]) +                                  
                                   sum(λ[m,t]*slackk[m,t] + penalty*v_p[m,t] for m = Mi, t = Tp) - 
                                   sum(ta[i]*stard1[i] for i = I) -
                                   sum(tb[i]*stard2[i,j,1] + tc[i]*stard2[i,j,2] for i = I, j = J[i]) -
                                   sum(λ[m,t]*sslackk[m,t] + penalty*sv_p[m,t] for m = Mi, t = Tp)
                )

    @objective(model, Min, total_tard)
   
    return model, slackk, v_p, b1, b2
end

function create_solve!(::Subproblem, d::JobShopProblem, I, λ)
    m, s, vp, b1, b2 = create_problem(Subproblem(), d, I, λ)
    optimize!(m)
    @show termination_status(m)
    @show objective_value(m)
    # finalize for CPLEX or EXPRESS
    GC.gc()
    return m, s, vp, b1, b2
end
