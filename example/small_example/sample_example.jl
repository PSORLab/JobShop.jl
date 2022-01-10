using Pkg, CPLEX

Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\Uploaded Jobshop\\Jobshop.jl")

using JobShop

jsprob = load_from_csv(@__DIR__)
jsprob.parameter.start_upper_bound = 8000.0
jsprob.parameter.ShiftLength = 18
jsprob.parameter.prob = 0.05
jsprob.parameter.prob_r = 0.2
jsprob.parameter.optimizer = CPLEX.Optimizer


sequential_solve!(jsprob)