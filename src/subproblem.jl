"""
$TYPEDSIGNATURES

A subproblem with `i ∈ I`.
"""
function create_problem(::Subproblem, jsprob::JobShopProblem, I::Vector{Int}, λ)
    @unpack J, T, R, p, ps, pr, w, d, U, O, M = jsprob
    model, b1, c1, b2, c2, bI1, bI2 = create_problem(SharedProblem(), jsprob, I)

    o = original_objective_ex!(m, jsprob, c1, c2, I)
    @variable(model, s[m ∈ M, t ∈ T] ≥ 0)
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
    @expression(model, ex_γ[i ∈ I, j ∈ J[i], m ∈ M, t ∈ T], sum(b1[i,j,k]   for i ∈ (t-p[i,j,m,1]):t))
    @expression(model, ex_τ[i ∈ I, j ∈ J[i], m ∈ M, t ∈ T], sum(b2[i,j,k,2] for i ∈ (t-p[i,j,m,1]):t))
    @expression(model, g[i ∈ I, m ∈ M, t ∈ T], sum(γ[i,j]*ex_γ[i,j,m,t] + τ[i,j]*ex_τ[i,j,m,t] for (k,j) in O[m] if k == i))
    @objective(model, Min, sum(o[i] + sum(λ[i,m,t]*(g[i,m,t] + s[m,t] - M[m]) for m ∈ M, t ∈ T) for i ∈ I))
    return model, o, g, s
end

"""
$TYPEDSIGNATURES

Updates a subproblem corresponding to (m,I) with new dual values.
"""
function update_problem!(::Subproblem, jsprob::JobShopProblem, i, λ)
    o, g, s, M, I, T, = jsprob.o[i], jsprob.g[i], jsprob.s[i], jsprob.M, jsprob.I[i], jsprob.T
    @objective(jsprob.m[i], Min, sum(o[i] + sum(λ[i,m,t]*(g[i,m,t] + s[m,t] - M[m,t]) for m ∈ M, t ∈ T) for i ∈ I))
    return
end
