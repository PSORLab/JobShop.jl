using Pkg
#=
Pkg.add(url="PATH_TO_JOBSHOP_PACKAGE HERE")
# or
=#
Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\JobShop.jl")

using JobShop

jsprob = load_from_csv(@__DIR__)
jsprob.initial_estimate = 1400.0
#solve!(jsprob)