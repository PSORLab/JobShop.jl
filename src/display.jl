function display_iteration(d::JobShopProblem, j)
    @unpack current_iteration, solve_time, heurestic_time = d.status
    if upper_bound[current_iteration] < 100000.0
        outstring = ""
        outstring *= @sprintf("Iteration:  %d, (S):  %d  ", current_iteration, i)
        outstring *= @sprintf("SDV:  %.3e  FC:  %.3e  ", lower_bound[current_iteration], upper_bound[current_iteration])
        outstring *= @sprintf("GAP:  %.3e  Qual.:  %.3e  ", current_abs_tol(d), Math.round(quality*1000)/1000)
        outstring *= @sprintf("Total Time: %.2e  %.2e  %.2e", Math.round(-(begin-end))/1000, solve_time, heurestic_time)
        outstring *= @sprintf("Norm:  %.3e  %.3e  %.3e  %.3e  %.3e	",norm, maxest, est, step, penalty)
    end
end