function display_iteration(d::JobShopProblem, i)
    @unpack maxest, current_estimate, penalty, current_iteration, start_time, solve_time, 
            heurestic_time, current_norm, current_step, lower_bound, upper_bound = d.status
    total_time = time() - d.status.start_time
    last_iteration_lbd = maximum(keys(lower_bound))
    last_iteration_ubd = maximum(keys(upper_bound))
    @show lower_bound
    @show upper_bound
    if (upper_bound[last_iteration_ubd] < 100000.0) || d.parameter.verbosity > 1
        outstring = ""
        outstring *= @sprintf("Iteration:  %d, (S):  %d  ", current_iteration, i)
        outstring *= @sprintf("SDV:  %.3e  FC:  %.3e  ", lower_bound[last_iteration_lbd], upper_bound[last_iteration_ubd])
        outstring *= @sprintf("GAP:  %.3e  Qual.:  %.3e  ", current_abs_gap(d), current_rel_gap(d))
        outstring *= @sprintf("Total Time: %.2e  %.2e  %.2e ", total_time, solve_time, heurestic_time)
        outstring *= @sprintf("Norm:  %.3e  %.3e  %.3e  %.3e  %.3e	", current_norm, maxest, current_estimate, current_step, penalty)
        println(outstring)
    end
end