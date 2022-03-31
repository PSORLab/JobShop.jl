
pos_sqr(x) = max(value(x), 0.0)^2

"""
$(TYPEDSIGNATURES)

Solve Lagrangian subproblem.
"""
function solve_subproblem(jsp::JobShopProblem, Ii::Vector{Int}, q::Int)

    subproblem_start = time()

    @unpack I, J, Jop, PartDue, MachineCap, MachineType, R, T, IJT, MIJ, 
            mult, sTard1, sTard2, sbI1, sbI2, sslackk, sv_p, sbI1_indx, 
            sbI2_indx, sbI1_indx_n, sbI2_indx_n, sbTime1, sbTime2, Tmax = jsp
    @unpack prob, prob_r, ShiftLength = jsp.parameter
    @unpack current_norm, current_iteration, prior_norm, prior_step, 
            current_step, lower_bound, estimate, penalty = jsp.status

    stime = time()

    m = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    configure!(Subproblem(), jsp, m)
    set_silent(m)

    @variable(m, 0 <= bTime1[i=Ii, j=Jop[i]] <= Tmax, Int, start = round(sbTime1[i,j]))
    @variable(m, 0 <= bTime2[i=Ii, j=Jop[i], j1=Jop[i], r=R] <= Tmax, Int, start=round(sbTime2[i,j,j1,r]))
    @variables(m, begin              
        -16 <= slackk[MachineType, T] <= 16
        0 <= slackk1[MachineType, T] <= 16
        0 <= v_p[MachineType, T] <= 16
        0 <= v_m[MachineType, T] <= 16
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

        0 <= y[i=Ii, Jop[i]] <= Tmax, Int
        0 <= cTime1[i=Ii, Jop[i]] <= Tmax, Int                 
        0 <= cTime2[i=Ii, Jop[i], Jop[i], R] <= Tmax, Int 
        0 <= ComTime1[Ii] <= Tmax, Int
        0 <= ComTime2[i=Ii, Jop[i], R] <= Tmax, Int
        bTimeI1[i=Ii, Jop[i], T], Bin
        bTimeI2[i=Ii, Jop[i], Jop[i], R, T], Bin                    
    end)
    
    @expressions(m, begin
        Tard1[i=I], W31[i]*2*T[end]
        Tard2[i=I, j=Jop[i], r=R], W32[i,j,r]*2*T[end]
        LB, sum(Tard1[i]*(1-prob)^J[i] for i in I) +
            sum(Tard2[i,j,1]*((1-prob)^(j-1) - (1-prob)^j)*(1-prob_r) for i in I, j in Jop[i]) + 
            sum(Tard2[i,j,2]*((1-prob)^(j-1) - (1-prob)^j)*prob_r for i in I, j in Jop[i]) +
            sum(mult[mt,t]*slackk[mt, t] for mt in MachineType, t in T) + 
            sum((penalty/100)*v_p[mt,t] for mt in MachineType, t in T)
        TotalTard,  sum(Tard1[i]*(1-prob)^J[i] for i in I) + 
                    sum(Tard2[i,j,1]*((1-prob)^(j-1) - (1-prob)^j)*(1-prob_r) for i in I, j in Jop[i]) + 
                    sum(Tard2[i,j,2]*((1-prob)^(j-1) - (1-prob)^j)*prob_r for i in I, j in Jop[i]) +
                    sum(mult[mt,t]*slackk[mt, t] for mt in MachineType, t in T) + 
                    sum((penalty/100)*v_p[mt,t] for mt in MachineType, t in T) -
                    (sum(sTard1[i]*(1-prob)^J[i] for i in I) +
                    sum(sTard2[i,j,1]*((1-prob)^(j-1) - (1-prob)^j)*(1-prob_r) for i in I, j in Jop[i]) + 
                    sum(sTard2[i,j,2]*((1-prob)^(j-1) - (1-prob)^j)*prob_r for i in I, j in Jop[i]) +
                    sum(mult[mt,t]*sslackk[mt, t] for mt in MachineType, t in T) + 
                    sum((penalty/100)*sv_p[mt,t] for mt in MachineType, t in T))
    end)

    @objective(m, Min, TotalTard)

    @constraints(m, begin
        [i=Ii, j=Jop[i]],   cTime1[i,j] + 1 <= bTime2[i,1,j,1]
        [i=Ii, j=Jop[i]],   cTime1[i,j] + 1 <= bTime2[i,j,j,2]

        W11 .<= alpha11 
        W21 .<= alpha21
        W31 .<= alpha31
        W11 .+ W21 .+ W31 .== 1 
        alpha11 .+ alpha31 .<= 1 
        alpha21 .== 1
        [i=Ii],   ComTime1[i] - PartDue[i] == W11[i]*(-T[end]) + W31[i]*2*T[end]

        [i=Ii, j=Jop[i], r=R],  W12[i,j,r] <= alpha12[i,j,r]
        [i=Ii, j=Jop[i], r=R],  W22[i,j,r] <= alpha22[i,j,r]
        [i=Ii, j=Jop[i], r=R],  W32[i,j,r] <= alpha32[i,j,r]
        [i=Ii, j=Jop[i], r=R],  W12[i,j,r] + W22[i,j,r] + W32[i,j,r] == 1
        [i=Ii, j=Jop[i], r=R],  alpha12[i,j,r] + alpha32[i,j,r] <= 1
        [i=Ii, j=Jop[i], r=R],  alpha22[i,j,r] == 1
        [i=Ii, j=Jop[i], r=R],  ComTime2[i,j,r] - PartDue[i] == W12[i,j,r]*(-T[end]) +  W32[i,j,r]*2*T[end]
    end)

    for i = Ii
        if J[i] >= 2
            @constraints(m, begin
                [j=1:(J[i]-1)],                       cTime1[i,j] + 1 <= bTime1[i,j+1]
                [j1=Jop[i], r=R, j=1:(J[i]-1)],  cTime2[i,j,j1,r] + 1 <= bTime2[i,j+1,j1,r]

                [j1=Jop[i]], y[i,j1] >= cTime1[i,j1]/ShiftLength
                [j1=Jop[i]], y[i,j1] <= cTime1[i,j1]/ShiftLength+0.9999
                [j1=Jop[i]], y[i,j1] <= (bTime2[i,1,j1,1]-1)/ShiftLength
                [j1=Jop[i]], y[i,j1] <= (bTime2[i,j1,j1,2]-1)/ShiftLength
            end)
        end
    end

    @constraints(m, begin
        [k=setdiff(I,Ii)],                    Tard1[k] == sTard1[k]
        [k=setdiff(I,Ii), j=Jop[k], r=R], Tard2[k,j,r] == sTard2[k,j,r]

        slackk .<= v_p
        slackk .>= -v_p
    end)
 
    s = zeros(length(MachineType), length(T))
    sbI1_i = sbI1_indx[q]
    sbI2_i = sbI2_indx[q]
    sbI1_in = sbI1_indx_n[q]
    sbI2_in = sbI2_indx_n[q]
    @inbounds for mi in MachineType, t in T
        temp = 0.0
        for (i,j,k) in sbI1_i[mi,t]
            temp += ((1-prob)^(j-1))*sbI1[i,j,k]
        end
        for (i,j,j1,k) in sbI2_i[mi,t]
            temp += ((1-prob)^(j-1)-(1-prob)^j)*(prob_r*sbI2[i,j,j1,2,k] + (1-prob_r)*sbI2[i,j,j1,1,k])
        end
        s[mi,t] = temp
    end

    @expression(m, mc_s1[mi=MachineType, t=T], sum(((1-prob)^(j-1))*bTimeI1[i,j,k] for (i,j,k) in sbI1_in[mi,t]))
    @expression(m, mc_s2[mi=MachineType, t=T], sum(((1-prob)^(j-1)-(1-prob)^(j))*(1-prob_r)*bTimeI2[i,j,j1,1,k] for (i,j,j1,k) in sbI2_in[mi,t]))
    @expression(m, mc_s3[mi=MachineType, t=T], sum(((1-prob)^(j-1)-(1-prob)^(j))*prob_r*bTimeI2[i,j,j1,2,k] for (i,j,j1,k) in sbI2_in[mi,t]))
    @constraint(m, [mi=MachineType, t=T], s[mi,t] + mc_s1[mi,t] + mc_s2[mi,t] + mc_s3[mi,t] - MachineCap[mi] + slackk1[mi,t] - slackk[mi,t] == 0.0) 

    @constraints(m, begin

        [i=Ii, j=Jop[i]], sum(sum(bTimeI1[i,j,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) == 1
        [i=Ii, j=Jop[i]], sum(sum(t*bTimeI1[i,j,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) == bTime1[i,j]
        [i=Ii, j=Jop[i]], sum(sum((t+Pro.t)*bTimeI1[i,j,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) - 1 == cTime1[i,j]
        
        [i=Ii, j=Jop[i], j1=Jop[i], r=R], sum(sum(bTimeI2[i,j,j1,r,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) == 1
        [i=Ii, j=Jop[i], j1=Jop[i], r=R], sum(sum(t*bTimeI2[i,j,j1,r,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) == bTime2[i,j,j1,r]
        [i=Ii, j=Jop[i], j1=Jop[i], r=R], sum(sum((t+Pro.t)*bTimeI2[i,j,j1,r,t] for t in T) for Pro in IJT if Pro.i == i && Pro.j == j) - 1 == cTime2[i,j,j1,r]
        
        [i=Ii],                 ComTime1[i]      == cTime1[i,J[i]]
        [i=Ii, j1=Jop[i], r=R], ComTime2[i,j1,r] == cTime2[i,J[i],j1,r]
    end)
    optimize!(m)

    jsp.status.time_solve_subprob += solve_time(m)
    valid_flag = valid_solve(Subproblem(), m)
    if valid_flag
        jsp.status.current_norm = sum(pos_sqr, slackk)
        jsp.status.lower_bound[current_iteration] = value(LB)
        jsp.status.lower_bound_time[current_iteration] = time() - jsp.status.time_start

        for mi in MachineType, t in T
            jsp.sslackk[mi,t] = value(slackk[mi,t])
            jsp.sv_p[mi,t] = value(v_p[mi,t])
        end
        jsp.mult .= max.(jsp.mult, 0.0)
        for t in T, i in Ii, j in Jop[i]
            jsp.sbI1[i,j,t] = value(bTimeI1[i,j,t])
            for j1 in Jop[i], r in R
                jsp.sbI2[i,j,j1,r,t] = value(bTimeI2[i,j,j1,r,t])
            end
        end
        for i in Ii
            jsp.sTard1[i] = value(Tard1[i])
            for j in Jop[i]
                jsp.sbTime1[i,j] = value(bTime1[i,j])
                for r in R
                    jsp.sTard2[i,j,r] = value(Tard2[i,j,r])
                    for j1 in Jop[i]
                        jsp.sbTime2[i,j,j1,r] = value(bTime2[i,j,j1,r])
                    end
                end
            end
        end
    end

    close_problem!(m)
    jsp.status.time_total_subprob = time() - subproblem_start
    return valid_flag
end
