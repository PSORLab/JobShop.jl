function solve_subproblem(jsp::JobShopProblem, Ii::Vector{Int})
    @unpack I, J, Jop, PartDue, MachineCap, MachineType, 
            R, T, IJT, MIJ, prob, prob_r, ShiftLength, mult, 
            sTard1, sTard2, sbI1, sbI2, sslackk, sv_p = jsp
    @unpack alpha_step = jsp.parameter
    @unpack current_norm, current_iteration, current_step, lower_bound, estimate, penalty = jsp.status

    m = direct_model(optimizer_with_attributes(jsp.parameter.optimizer))
    configure!(Subproblem(), jsp, m)
    set_silent(m)

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

        0 <= y[i=Ii, Jop[i]] <= 1000, Int
        0 <= cTime1[i=Ii, Jop[i]] <= 1000, Int
        0 <= bTime1[i=Ii, Jop[i]] <= 1000, Int                   
        0 <= cTime2[i=Ii, Jop[i], Jop[i], R] <= 1000, Int 
        0 <= bTime2[i=Ii, Jop[i], Jop[i], R] <= 1000, Int
        0 <= ComTime1[Ii] <= 1000, Int
        0 <= ComTime2[i=Ii, Jop[i], R] <= 1000, Int
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

    @constraints(m, begin
        [i in Ii],                               ComTime1[i] - cTime1[i,J[i]] == 0.0
        [i=Ii, j1=Jop[i], r=R],                  ComTime2[i,j1,r] == cTime2[i,J[i],j1,r]
        slackk .<= v_p
        slackk .>= -v_p
        [i=Ii, j=Jop[i]],                        cTime1[i,j] + 1 <= bTime2[i,1,j,1]
        [i=Ii, j=Jop[i]],                        cTime1[i,j] + 1 <= bTime2[i,j,j,2]
        [i=Ii, j=1:(J[i]-1)],                    cTime1[i,j] + 1 <= bTime1[i,j+1]
        [i=Ii],                                  ComTime1[i]     == W11[i]*(-T[end]) + W31[i]*2*T[end] + PartDue[i]
        [i=Ii],                                  ComTime1[i]     == cTime1[i,J[i]]
        [i=Ii, j=Jop[i],  r=R],          ComTime2[i,j,r] - PartDue[i]  == W12[i,j,r]*(-T[end]) + W32[i,j,r]*2*T[end]
        [i=Ii, j1=Jop[i], r=R, j = 1:(J[i]-1)],  cTime2[i,j,j1,r] <= bTime2[i,j+1,j1,r] - 1
        [k ∈ setdiff(I,Ii)],                            Tard1[k] == sTard1[k]
        [k ∈ setdiff(I,Ii), j = Jop[k], r = R], Tard2[k,j,r] == sTard2[k,j,r]
        [i = Ii, j = Jop[i]], sum(bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i = Ii, j = Jop[i]], sum(t*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime1[i,j]
        [i = Ii, j = Jop[i]], sum((t+P.t)*bTimeI1[i,j,t] for t in T, P in IJT if has_ij(P,i,j)) - 1 == cTime1[i,j]
        W11 .<= alpha11 
        W21 .<= alpha21
        W31 .<= alpha31
        W11 .+ W21 .+ W31 .== 1 
        alpha11 .+ alpha31 .<= 1 
        alpha21 .== 1
        [i=I, j=Jop[i], r=R], W12[i,j,r] <= alpha12[i,j,r]
        [i=I, j=Jop[i], r=R], W22[i,j,r] <= alpha22[i,j,r]
        [i=I, j=Jop[i], r=R], W32[i,j,r] <= alpha32[i,j,r]
        [i=I, j=Jop[i], r=R], W12[i,j,r] + W22[i,j,r] + W32[i,j,r] == 1
        [i=I, j=Jop[i], r=R], alpha12[i,j,r] + alpha32[i,j,r] <= 1
        [i=I, j=Jop[i], r=R], alpha22[i,j,r] == 1
        [i = Ii, j = Jop[i], j1 = Jop[i], r = R], sum(bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == 1
        [i = Ii, j = Jop[i], j1 = Jop[i], r = R], sum(t*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) == bTime2[i,j,j1,r]
        [i = Ii, j = Jop[i], j1 = Jop[i], r = R], sum((t+P.t)*bTimeI2[i,j,j1,r,t] for t in T, P in IJT if has_ij(P,i,j)) - 1  == cTime2[i,j,j1,r]
    end)

    for i = Ii
        if J[i] >= 2
            @constraints(m, begin
                [j1=Jop[i]], y[i,j1] >= cTime1[i,j1]/ShiftLength
                [j1=Jop[i]], y[i,j1] <= cTime1[i,j1]/ShiftLength+0.9999
                [j1=Jop[i]], y[i,j1] <= (bTime2[i,1,j1,1]-1)/ShiftLength
                [j1=Jop[i]], y[i,j1] <= (bTime2[i,j1,j1,2]-1)/ShiftLength
            end)
        end
    end

    # fast unpack for surrogate values
    outi = fill(true,length(I))
    for i in Ii
        outi[i] = false
    end
    s1 = zeros(length(MachineType),length(T))
    for mi in MachineType, t in T
        t2 = 0.0
        for Ma in MIJ
            c = (1-prob)^(Ma.j-1)
            t1 = 0.0
            for p in IJT
                if (p.i == Ma.i) && (p.j == Ma.j) && (Ma.m == mi) && outi[p.i] # replace this !(p.part in Ii)
                    for k in (t - p.t + 1):t
                        if t - p.t + 1 >= 1
                            t1 += sbI1[Ma.i,Ma.j,k]
                        end
                    end
                end
            end
            t2 += c*t1
        end
        s1[mi,t] = t2
    end
    s2 = zeros(length(MachineType),length(T))
    s3 = zeros(length(MachineType),length(T))
    for mi in MachineType, t in T
        t2 = 0.0
        t3 = 0.0
        for Ma in jsp.MIJ, p in jsp.IJT, j in Jop[Ma.i]
            if (p.i == Ma.i) && (p.j == Ma.j) && (Ma.m == mi) && !(p.i in Ii)
                t2a = 0.0
                t3a = 0.0
                for k in (t - p.t+1):t
                    if t - p.t+1 >= 1
                        t2a += sbI2[Ma.i,Ma.j,j,1,k]
                        t3a += sbI2[Ma.i,Ma.j,j,2,k]
                    end
                end
                t2 += ((1-prob)^(j-1)-(1-prob)^(j))*(1-prob_r)*t2a
                t3 += ((1-prob)^(j-1)-(1-prob)^(j))*prob_r*t3a
            end
        end
        #@show mi, t, indx_c
        s2[mi,t] = t2
        s3[mi,t] = t3
    end
    @constraint(m, [mi in MachineType, t in T],
        s1[mi,t] + s2[mi,t] + s3[mi,t] +
        sum(((1-prob)^(Ma.j-1))*sum(bTimeI1[Ma.i,Ma.j,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for p in IJT, Ma in MIJ if (p.i == Ma.i && p.j == Ma.j && Ma.m == mi && (p.i in Ii))) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*(1-prob_r)*sum(bTimeI2[Ma.i,Ma.j,j,1,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if p.i == Ma.i && p.j == Ma.j && Ma.m == mi && (p.i in Ii)) +
        sum(((1-prob)^(j-1)-(1-prob)^(j))*prob_r*sum(bTimeI2[Ma.i,Ma.j,j,2,k] for k in (t - p.t+1):t if (t - p.t+1 >= 1)) for Ma in MIJ, p in IJT, j in Jop[Ma.i] if p.i == Ma.i && p.j == Ma.j && Ma.m == mi && (p.i in Ii))
         - MachineCap[mi] + slackk1[mi,t] == slackk[mi,t])
    
    @objective(m, Min, TotalTard)
    
    optimize!(m)

    jsp.status.solve_time += solve_time(m)
    valid_flag = valid_solve(Subproblem(), m)
    if valid_flag
        jsp.status.current_norm = sum(x -> max(value(x), 0.0)^2, slackk)
        jsp.status.lower_bound[current_iteration] = value(LB)
        if (current_iteration > 28) && (estimate < 100000) && (estimate - value(LB) > 0)
            jsp.status.current_step = alpha_step/10.0*(estimate - value(LB))/jsp.status.current_norm
        end
        for mi in MachineType, t in T
            temp = value(slackk[mi,t])
            jsp.mult[mi,t] += current_step*temp
            jsp.sslackk[mi,t] = temp
            jsp.sv_p[mi,t] = value(v_p[mi,t])
        end
        jsp.mult .= max.(jsp.mult, 0.0)
        for t in T, i in Ii, j in Jop[i]
            jsp.sbI1[i,j,t] = value(bTimeI1[i,j,t])
            for j1 in Jop[i], r in R
                jsp.sbI2[i,j,j1,r,t] = value(bTimeI2[i,j,j1,r,t])
            end
        end
        for i in Ii, j in Jop[i], r in R
            jsp.sTard1[i]     = value(Tard1[i])
            jsp.sTard2[i,j,r] = value(Tard2[i,j,r])
            jsp.sbTime1[i,j]  = value(bTime1[i,j])
        end
    end
    return valid_flag
end
