"""
$TYPEDSIGNATURES

Constructs the feasibility problem called after convergence using 
estimates for  starting times obtain by Lagrangian relaxation.
"""
function solve_problem(::FeasibilityProblem, jsp::JobShopProblem)

    feasibility_start = time()

    @unpack I, J, Jop, PartDue, MachineCap, MachineType, R, T, IJT, MIJ, sbTime1 = jsp
    @unpack upper_bound = jsp.status
    @unpack prob, prob_r, ShiftLength, feasibility_window = jsp.parameter 

    m = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    configure!(FeasibilityProblem(), jsp, m)
    set_silent(m)

    @variables(m, begin          
        0 <= cTime1[i=I, Jop[i]] <= 1000, Int
        0 <= cTime2[i=I, Jop[i], Jop[i], R] <= 1000, Int 
        bTimeI1[i=I, Jop[i], T], Bin

        0 <= bTime1[i=I, Jop[i]] <= 1000, Int                   
        0 <= bTime2[i=I, Jop[i], Jop[i], R] <= 1000, Int
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
        [i=I, j1=Jop[i], r=R, j = 1:(J[i]-1)],  cTime2[i,j,j1,r] <= bTime2[i,j+1,j1,r] - 1  # TODO: maybe...

        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] >= cTime1[i,j1]/ShiftLength
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= cTime1[i,j1]/ShiftLength+0.9999
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= (bTime2[i,1,j1,1]-1)/ShiftLength
        [i=I, j1=Jop[i]; J[i] >= 2],    y[i,j1] <= (bTime2[i,j1,j1,2]-1)/ShiftLength
    end)

    @constraint(m, [mi in MachineType, t in T],
        sum(((1-prob)^(Ma.j-1))*sum(bTimeI1[Ma.i,Ma.j,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for p in IJT, Ma in MIJ if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi)) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*(1-prob_r)*sum(bTimeI2[Ma.i,Ma.j,j,1,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi)) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*prob_r*sum(bTimeI2[Ma.i,Ma.j,j,2,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi))
        <= MachineCap[mi] # TODO: maybe...
    )

    @constraints(m, begin
        [i=I, j=Jop[i]], bTime1[i,j] - sbTime1[i,j] <= feasibility_window
        [i=I, j=Jop[i]], bTime1[i,j] - sbTime1[i,j] >= -feasibility_window

        [i=I, j=Jop[i]], sum(bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i=I, j=Jop[i]], sum(t*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime1[i,j]
        [i=I, j=Jop[i]], sum((t+P.t)*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) - 1 == cTime1[i,j]

        [i=I, j=Jop[i], j1=Jop[i], r=R], sum(bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i=I, j=Jop[i], j1=Jop[i], r=R], sum(t*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime2[i,j,j1,r]
        [i=I, j=Jop[i], j1=Jop[i], r=R], sum((t+P.t)*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) - 1  == cTime2[i,j,j1,r]

        [i=I, j1=Jop[i], r=R],  ComTime2[i,j1,r] == cTime2[i,J[i],j1,r]
        [i=I],                  ComTime1[i]      == cTime1[i,J[i]]
    end)

    optimize!(m)
    @show termination_status(m)
    jsp.status.time_solve_feasibility += solve_time(m)
    valid_flag = valid_solve(FeasibilityProblem(), m)
    if valid_flag
        @show objective_value(m)
        jsp.status.upper_bound[jsp.status.current_iteration] = objective_value(m)
    end

    close_problem!(m)
    jsp.status.time_total_feasibility = time() - feasibility_start

    return valid_flag
end

"""
$TYPEDSIGNATURES

Checks to see whether the feasibility problem should be solved.
"""
function use_problem(::FeasibilityProblem, d::JobShopProblem)
    @unpack current_norm, current_iteration = d.status
    @unpack feasible_norm_limit, feasible_interval, feasible_start, verbosity = d.parameter
    flag = current_iteration >= feasible_start
    flag &= iszero(mod(current_iteration, feasible_interval))
    flag &= current_norm < feasible_norm_limit
    (verbosity > 1) && (flag && println("Feasibility problem will be solved."))
    return flag
end