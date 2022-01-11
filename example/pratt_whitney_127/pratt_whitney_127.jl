using Pkg, CPLEX, JLD
#=
Pkg.add(url="PATH_TO_JOBSHOP_PACKAGE HERE")
# or
=#
#Pkg.rm("JobShop")
Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\Uploaded Jobshop\\Jobshop.jl")

using JobShop

jsprob = load_from_csv(@__DIR__)
jsprob.parameter.start_upper_bound = 8000.0
jsprob.parameter.ShiftLength = 18
jsprob.parameter.prob = 0.05
jsprob.parameter.prob_r = 0.2
jsprob.T = 1:220
jsprob.parameter.start_norm = 100.0
jsprob.parameter.start_step = 0.006
jsprob.parameter.penalty = 120.0
jsprob.parameter.optimizer = CPLEX.Optimizer
jsprob.parameter.alpha_step_1 = 0.5/20
jsprob.parameter.alpha_step_2 = 0.5
jsprob.parameter.feasible_interval = 7
jsprob.parameter.feasible_norm_limit = 3.0
jsprob.parameter.stepsize_interval = 25
jsprob.parameter.stepsize_start = 37.5

sequential_solve!(jsprob)
save(joinpath(@__DIR__, "127part.jld", "jsp", jsprob))