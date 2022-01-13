using Pkg, CPLEX, JLD

Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\Uploaded Jobshop\\Jobshop.jl")

using JobShop

# define problem
jsprob = load_from_csv(@__DIR__)
jsprob.T = 1:75
jsprob.parameter.ShiftLength = 18
jsprob.parameter.prob = 0.05
jsprob.parameter.prob_r = 0.2

# set hyperparameter
jsprob.parameter.start_upper_bound = 8000.0
jsprob.parameter.start_norm = 100.0
jsprob.parameter.start_step = 0.006
jsprob.parameter.optimizer = CPLEX.Optimizer
jsprob.parameter.penalty = 60.0
jsprob.parameter.alpha_step_1 = 0.5/20
jsprob.parameter.alpha_step_2 = 0.5
jsprob.parameter.feasible_interval = 21
jsprob.parameter.feasible_norm_limit = 0.7 # TODO: CHANGE BACK TO 2
jsprob.parameter.feasibility_window = 2
jsprob.parameter.stepsize_interval = 40
jsprob.parameter.stepsize_start = 37.5
jsprob.parameter.iteration_limit = 1000

# solve
sequential_solve!(jsprob)
#save(joinpath(@__DIR__, "simple_example.jld", "jsp", jsprob))