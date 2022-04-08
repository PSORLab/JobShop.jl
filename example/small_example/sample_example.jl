using Pkg, CPLEX, JuMP, JLD

#Pkg.develop(path="C:\\Users\\maw16110\\Desktop\\Jobshop.jl\\Jobshop.jl")

using JobShop

# define problem
jsprob = load_from_csv(@__DIR__)
jsprob.T = 1:75
jsprob.Tmax = 1000
jsprob.parameter.ShiftLength = 18
jsprob.parameter.prob = 0.05
jsprob.parameter.prob_r = 0.2

# set hyperparameter
jsprob.parameter.random_seed = 0
jsprob.parameter.start_upper_bound = 8000.0
jsprob.parameter.start_norm = 125.0
jsprob.parameter.start_step = 0.05
jsprob.parameter.optimizer = CPLEX.Optimizer
jsprob.parameter.penalty = 60.0
jsprob.parameter.feasible_norm_limit = 2.5
jsprob.parameter.feasible_solve_count = 4
jsprob.parameter.feasible_start = 50
jsprob.parameter.iteration_limit = 10000000000
jsprob.parameter.penalty_iteration = 10000000
jsprob.parameter.penalty_factor = 1.0
jsprob.parameter.feasible_solve_time = 120.0
jsprob.parameter.verbosity = 1
jsprob.parameter.use_stepsize_program = true

# add custom starting time restriction on feasibility problem
function feasibility_starting_time_constraints!(::Ext, jsp::JobShopProblem, m, bTime1, sbTime1)
    I = jsp.I
    @constraints(m, begin
        [i=I, j=1:3], bTime1[i,j] - sbTime1[i,j] <= 3
        [i=I, j=1:3], bTime1[i,j] - sbTime1[i,j] >= -3
    end)
    return
end

# solve
sequential_solve!(jsprob)

 # save result
save(joinpath(@__DIR__, "simple_example_new.jld"), "jsp", jsprob)
