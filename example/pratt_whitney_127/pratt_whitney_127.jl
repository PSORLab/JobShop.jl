using Pkg, CPLEX, JLD
#=
Pkg.add(url="PATH_TO_JOBSHOP_PACKAGE HERE")
# or
=#
#Pkg.rm("JobShop")
#Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\Jobshop.jl\\Jobshop.jl")

using JobShop

function solve_example(pi, pf, sc, fw, nl)
    jsprob = load_from_csv(@__DIR__)
    jsprob.parameter.start_upper_bound = 8000.0
    jsprob.parameter.ShiftLength = 18
    jsprob.parameter.prob = 0.05
    jsprob.parameter.prob_r = 0.2
    jsprob.T = 1:220

    jsprob.parameter.start_norm = 200.0
    jsprob.parameter.start_step = 0.006
    jsprob.parameter.penalty = 120.0
    jsprob.parameter.optimizer = CPLEX.Optimizer

    jsprob.parameter.feasible_norm_limit = nl
    jsprob.parameter.feasibility_window = fw
    jsprob.parameter.feasible_solve_count = sc
    jsprob.parameter.iteration_limit = 10000
    jsprob.parameter.penalty_iteration = pi
    jsprob.parameter.penalty_factor = pf
    jsprob.parameter.feasible_solve_time = 120.0
    jsprob.parameter.random_seed = 1

    sequential_solve!(jsprob)

    save(joinpath(@__DIR__, "127part.jld"), "jsp", jsprob)
end


solve_example(1, 1.0, 6, 3, 2.0)