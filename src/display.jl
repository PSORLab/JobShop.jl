function display_iteration(d::JobShopProblem, j)
    @unpack current_iteration, start_time, solve_time, heurestic_time, current_norm, current_step = d.status
    total_time = time() - d.start_time
    #maxest = 
    #est =  
    #penalty = 
    if upper_bound[current_iteration] < 100000.0
        outstring = ""
        outstring *= @sprintf("Iteration:  %d, (S):  %d  ", current_iteration, i)
        outstring *= @sprintf("SDV:  %.3e  FC:  %.3e  ", lower_bound[current_iteration], upper_bound[current_iteration])
        outstring *= @sprintf("GAP:  %.3e  Qual.:  %.3e  ", current_abs_gap(d), current_rel_gap(d))
        outstring *= @sprintf("Total Time: %.2e  %.2e  %.2e", total_time, solve_time, heurestic_time)
        outstring *= @sprintf("Norm:  %.3e  %.3e  %.3e  %.3e  %.3e	", current_norm, maxest, est, current_step, penalty)
        println(outstring)
    end
end