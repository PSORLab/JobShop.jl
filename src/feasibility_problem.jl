
function feasibility_starting_time_constraints!(ext::AbstractExt, jsp::JobShopProblem, m, bTime1, sbTime1)
    @unpack I, J, Jop = jsp
    @unpack feasibility_window = jsp.parameter 
    @constraints(m, begin
        [i=I, j=Jop[i]], bTime1[i,j] - sbTime1[i,j] <= feasibility_window + 0.001
        [i=I, j=Jop[i]], bTime1[i,j] - sbTime1[i,j] >= -feasibility_window - 0.001
    end)
    return
end

"""
$(TYPEDSIGNATURES)

Constructs the feasibility problem called after convergence using 
estimates for the starting times obtain by Lagrangian relaxation.
"""
function solve_problem(::FeasibilityProblem, jsp::JobShopProblem)

    feasibility_start = time()

    @unpack I, J, Jop, PartDue, MachineCap, MachineType, R, T, IJT, MIJ, 
            sbTime1, sbTime2, mult, sslackk, sv_p = jsp
    @unpack upper_bound, penalty = jsp.status
    @unpack prob, prob_r, ShiftLength, feasible_solve_count, feasible_solve_time, ext = jsp.parameter 

    m = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    configure!(FeasibilityProblem(), jsp, m)
    set_silent(m)
    set_time_limit_sec(m, feasible_solve_time)

    @variable(m, 0 <= bTime1[i=I, j=Jop[i]] <= 1000, Int, start = round(sbTime1[i,j]))
    @variable(m, 0 <= bTime2[i=I, j=Jop[i], j1=Jop[i], r=R] <= 1000, Int, start=round(sbTime2[i,j,j1,r]))
    @variables(m, begin          
        0 <= cTime1[i=I, Jop[i]] <= 1000, Int
        0 <= cTime2[i=I, Jop[i], Jop[i], R] <= 1000, Int
        bTimeI1[i=I, Jop[i], T], Bin
        bTimeI2[i=I, Jop[i], Jop[i], R, T], Bin

        0 <= ComTime1[I] <= 1000, Int
        0 <= ComTime2[i=I,Jop[i],R] <= 1000, Int

        0 <= y[i=I, Jop[i]] <= 1000, Int   

        # variables used for first SOS rearrangement
        0 <= W11[I] <= 1
        0 <= W21[I] <= 1
        0 <= W31[I] <= 1
        alpha11[I], Bin  
        alpha21[I], Bin 
        alpha31[I], Bin

        # variables used for second SOS rearrangement
        0 <= W12[i=I, j=Jop[i], R] <= 1
        0 <= W22[i=I, j=Jop[i], R] <= 1
        0 <= W32[i=I, j=Jop[i], R] <= 1
        alpha12[i=I, j=Jop[i], R], Bin
        alpha22[i=I, j=Jop[i], R], Bin
        alpha32[i=I, j=Jop[i], R], Bin
     end)

    @expressions(m, begin
        Tard1[i=I], W31[i]*2*T[end]
        Tard2[i=I, j=Jop[i], r=R], W32[i,j,r]*2*T[end]
    end)

    @expression(m, TotalTard, 
        sum(Tard1[i]*(1-prob)^J[i] for i in I) + 
        sum(Tard2[i,j,1]*((1-prob)^(j-1) - (1-prob)^j)*(1-prob_r) for i in I, j in Jop[i]) + 
        sum(Tard2[i,j,2]*((1-prob)^(j-1) - (1-prob)^j)*prob_r for i in I, j in Jop[i])
        )

    @objective(m, Min, TotalTard)

    @constraints(m, begin

        [i=I, j=Jop[i]],        cTime1[i,j] + 1 <= bTime2[i,1,j,1]
        [i=I, j=Jop[i]],        cTime1[i,j] + 1 <= bTime2[i,j,j,2]

        W11 .<= alpha11 
        W21 .<= alpha21
        W31 .<= alpha31
        W11 .+ W21 .+ W31 .== 1 
        alpha11 .+ alpha31 .<= 1 
        alpha21 .== 1
        [i=I],  ComTime1[i] - PartDue[i] == W11[i]*(-T[end]) + W31[i]*2*T[end]

        [i=I, j=Jop[i], r=R],   W12[i,j,r] <= alpha12[i,j,r]
        [i=I, j=Jop[i], r=R],   W22[i,j,r] <= alpha22[i,j,r]
        [i=I, j=Jop[i], r=R],   W32[i,j,r] <= alpha32[i,j,r]
        [i=I, j=Jop[i], r=R],   W12[i,j,r] + W22[i,j,r] + W32[i,j,r] == 1
        [i=I, j=Jop[i], r=R],   alpha12[i,j,r] + alpha32[i,j,r] <= 1
        [i=I, j=Jop[i], r=R],   alpha22[i,j,r] == 1
        [i=I, j=Jop[i], r=R],   ComTime2[i,j,r] - PartDue[i] == W12[i,j,r]*(-T[end]) + W32[i,j,r]*2*T[end]

        [i=I, j=1:(J[i]-1)],                    cTime1[i,j] + 1  <= bTime1[i,j+1]
        [i=I, j1=Jop[i], r=R, j = 1:(J[i]-1); J[i] >= 2],  cTime2[i,j,j1,r] + 1 <= bTime2[i,j+1,j1,r]

        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] >= cTime1[i,j1]/ShiftLength
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= cTime1[i,j1]/ShiftLength + 0.9999
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= (bTime2[i,1,j1,1]-1)/ShiftLength
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= (bTime2[i,j1,j1,2]-1)/ShiftLength
    end)

    @constraint(m, [mi in MachineType, t in T],
        sum(((1-prob)^(Ma.j-1))*sum(bTimeI1[Ma.i,Ma.j,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for p in IJT, Ma in MIJ if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi)) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*(1-prob_r)*sum(bTimeI2[Ma.i,Ma.j,j,1,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi)) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*prob_r*sum(bTimeI2[Ma.i,Ma.j,j,2,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi))
        <= MachineCap[mi] + 0.00001 # TODO: maybe...
    )

    @constraints(m, begin
        [i=I, j=Jop[i]], sum(bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i=I, j=Jop[i]], sum(t*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime1[i,j]
        [i=I, j=Jop[i]], sum((t+P.t)*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) - 1 == cTime1[i,j]

        [i=I, j=Jop[i], j1=Jop[i], r=R], sum(bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i=I, j=Jop[i], j1=Jop[i], r=R], sum(t*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime2[i,j,j1,r]
        [i=I, j=Jop[i], j1=Jop[i], r=R], sum((t+P.t)*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) - 1  == cTime2[i,j,j1,r]

        [i=I, j1=Jop[i], r=R],  ComTime2[i,j1,r] == cTime2[i,J[i],j1,r]
        [i=I],                  ComTime1[i]      == cTime1[i,J[i]]
    end)

    feasibility_starting_time_constraints!(ext, jsp, m, bTime1, sbTime1)

    valid_flag = false
    for i=1:feasible_solve_count
        optimize!(m)
        jsp.status.time_solve_feasibility += solve_time(m)
        valid_flag = valid_solve(FeasibilityProblem(), m)
        if valid_flag
            if primal_status(m) == MOI.FEASIBLE_POINT
                bc_i = jsp.status.current_iteration + i - 1
                bc_t = time() - jsp.status.time_start
                jsp.status.lower_bound_time[bc_i] = bc_t
                jsp.status.upper_bound_time[bc_i] = bc_t
                jsp.status.lower_bound[bc_i] = jsp.status.lower_bound[maximum(keys(jsp.status.lower_bound))]
                jsp.status.upper_bound[bc_i] = objective_value(m)
            end
        end
    end
    jsp.status.feasible_problem_found = valid_flag
    println("Was a feasible point found? $(valid_flag)")

    for i=I, j=Jop[i]
        jsp.feasible_bTime1[(i,j)] = round(value(bTime1[i,j]))
    end
    for i=I, j1 =Jop[i], j2 = Jop[i], r = R
        jsp.feasible_bTime2[(i,j1,j2,r)] = round(value(bTime2[i,j1,j2,r]))
    end

    close_problem!(m)
    jsp.status.time_total_feasibility = time() - feasibility_start

    return valid_flag
end

"""
$(TYPEDSIGNATURES)

Checks to see whether the feasibility problem should be solved.
"""
function use_problem(::FeasibilityProblem, d::JobShopProblem)
    @unpack current_norm, current_iteration = d.status
    @unpack feasible_norm_limit, feasible_start = d.parameter
    flag = current_iteration >= feasible_start
    flag &= current_norm < feasible_norm_limit
    return flag
end