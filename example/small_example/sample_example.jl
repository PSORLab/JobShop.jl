using Pkg, CPLEX, JuMP, JLD, UnPack

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
jsprob.parameter.start_upper_bound = 100000.0
jsprob.parameter.start_norm = 100.0
jsprob.parameter.start_step = 0.1
jsprob.parameter.optimizer = CPLEX.Optimizer
jsprob.parameter.penalty = 75.0
jsprob.parameter.feasible_norm_limit = 1.5
jsprob.parameter.feasible_solve_count = 5
jsprob.parameter.feasibility_window = 6
jsprob.parameter.feasible_start = 200
jsprob.parameter.iteration_limit = 10000000000
jsprob.parameter.penalty_iteration = 10000000
jsprob.parameter.penalty_factor = 1.05
jsprob.parameter.feasible_solve_time = 60.0
jsprob.parameter.verbosity = 1
jsprob.parameter.use_stepsize_program = true

# add custom starting time restriction on feasibility problem
function JobShop.feasibility_starting_time_constraints!(::Ext, jsp::JobShopProblem, m, bTime1, sbTime1)
    @unpack I, J, Jop = jsp
    @unpack feasibility_window = jsp.parameter
    lt_dict = Dict{Any,Int}()
    gt_dict = Dict{Any,Int}()
    for i in I
        lt_dict[i,1] = 1
        gt_dict[i,1] = 1
        lt_dict[i,2] = 1
        gt_dict[i,2] = 1
        lt_dict[i,3] = 3
        gt_dict[i,3] = 3
        lt_dict[i,4] = 100
        gt_dict[i,4] = 100
        lt_dict[i,5] = 100
        gt_dict[i,5] = 100
    end
    @constraints(m, begin
        fw_lt[i=I, j=1:5], bTime1[i,j] <= lt_dict[i,j] + sbTime1[i,j] + 0.001
        fw_gt[i=I, j=1:5], bTime1[i,j] >= -gt_dict[i,j] + sbTime1[i,j] - 0.001
    end)
    return fw_lt, fw_gt
end

# solve
sequential_solve!(jsprob)

 # save result
save(joinpath(@__DIR__, "simple_example_new.jld"), "jsp", jsprob)
