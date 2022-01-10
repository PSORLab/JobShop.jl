"""
$TYPEDSIGNATURES

Displays the statistics for each iteration.
"""
function display_iteration(j::JobShopProblem, i::Int)
    @unpack maxest, estimate, penalty, current_iteration, start_time, solve_time, 
            heurestic_time, current_norm, current_step = j.status
    total_time = time() - start_time
    if (current_upper_bound(j) < 100000.0) || (j.parameter.verbosity > 1)
        outstring = @sprintf("Iteration:  %d, (S):  %d  ", current_iteration, i)
        outstring *= @sprintf("SDV:  %.3e  FC:  %.3e  ", current_lower_bound(j), current_upper_bound(j))
        outstring *= @sprintf("GAP:  %.3e  Qual.:  %.3e  ", current_abs_gap(j), current_rel_gap(j))
        outstring *= @sprintf("Total Time: %.2e  %.2e  %.2e ", total_time, solve_time, heurestic_time)
        outstring *= @sprintf("Norm:  %.3e  %.3e  %.3e  %.3e  %.3e	", current_norm, maxest, estimate, current_step, penalty)
        println(outstring)
    end
end