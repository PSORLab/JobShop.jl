using Pkg, CPLEX
#=
Pkg.add(url="PATH_TO_JOBSHOP_PACKAGE HERE")
# or
=#
#Pkg.rm("JobShop")
Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\Uploaded Jobshop\\Jobshop.jl")

using JobShop

jsprob = load_from_csv(@__DIR__)
jsprob.parameter.start_estimate = 1400.0
jsprob.parameter.start_upper_bound = 8000.0
jsprob.parameter.optimizer = CPLEX.Optimizer


sequential_solve!(jsprob)