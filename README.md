# Jobshop.jl
Surrogate Lagrangian Jobshop Solver

## Problem Specification
The jobshop scheduling problem is specified by a `jsp::JobShopProblem` mutable structure. Currently, persistent 
storage and specification of the problem is done using a series of .csv files which are loaded into intermediate 
`DataFrames` objects then used to populate the `jsp` structure. In total eight .csv files are used which are 
loaded into dataframes then used:
- **part_due**: Contains a column `due`. The i-th row of column `due` has the due date of the i-th part.  
- **part_group**: Contains columns `part` and `group`. Defines the Lagrangian subproblem (group) that each part belongs to. 
- **part_operation_num**: Contains a column `num` which is the number of operations required to complete the i-th part which corresponds to the i-th row.
- **part_operation_rework**: Currently, unused. Will be used to set rework propababilities specific to each part and operation rather than a single rework probability for all part and operation combinations.
- **part_operation_scrap**: Currently, unused. Will be used to set scrap propababilities specific to each part and operation rather than a single scrap probability for all part and operation combinations.
- **part_operation_time**: Contains three columns `part`, `op`, and `time` used to specific the time taken to process each part on a given operation.
- **machine_part_op**: Describes the part and operation combinations that may be processes on a given machine.
- **machine_capacity**: Describes the capacity for each machine where the ith entry in the `capacity` column corresponds to the capacity of machine `i`.

The following parameters must also be set in order to fully the problem.

- **ShiftLength**: Typical shift length considered.
- **prob**: Set in the `parameter` subfield of `jsp`. Probability of scrap. 
- **prob_r**: Set in the `parameter` subfield of `jsp`. Probability of rework.
- **T**: UnitRange defining the number of timesteps considered.
- **Tmax**: Maximum time for scheduling problem.

## Solver Metaparameters
Options set as fields of the `parameter` subfield.

- **start_upper_bound**: Initial estimate of upper bounds.
- **start_norm**: Starting norm.
- **start_step**: Starting step size.
- **penalty**: Starting penalty.
- **optimizer**: Sets the optimizer used to solve subproblems and the feasibility problem.
- **feasible_norm_limit**: Below this absolute norm limit, the algorithm will start to try to solve the feasibility problem.
- **feasibility_window**: The neighbor about the Lagrangian subproblem solution a feasible solution will be searched for.
- **feasible_solve_count**: Number of feasibility problems to solve in sequence to generate feasible upper bounds.
- **iteration_limit**: Maximum number of subproblem iterations used.
- **penalty_iteration**: Number of iterations after which the penalty will be updated.
- **penalty_factor**: Factor used to update the penalty.
- **feasible_solve_time**: Solve time allowed each feasible solve iteration performed.
- **random_seed**: Random seed used for all routines called in the Random module of julia.

