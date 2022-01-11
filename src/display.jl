"""
$TYPEDSIGNATURES

Displays the statistics for each iteration.
"""
function display_iteration(j::JobShopProblem, i::Int)
    @unpack maxest, estimate, penalty, current_iteration, time_start, time_solve_subprob, 
            time_solve_stepsize, current_norm, current_step = j.status
    total_time = time() - time_start
    if (current_upper_bound(j) < 100000.0) || (j.parameter.verbosity > 1)
        outstring = @sprintf("Iteration:  %4i, (S):  %4i  ", current_iteration, i)
        outstring *= @sprintf("SDV:  %.3e  FC:  %.3e  ", current_lower_bound(j), current_upper_bound(j))
        outstring *= @sprintf("GAP:  %.3e  Qual.:  %.3e  ", current_abs_gap(j), current_rel_gap(j))
        outstring *= @sprintf("Total Time: %.1f  %.1f  %.1f ", total_time, time_solve_subprob, time_solve_stepsize)
        outstring *= @sprintf("Norm:  %.3e  %.3e  %.3e  %.3e  %.3e	", current_norm, maxest, estimate, current_step, penalty)
        println(outstring)
    end
end