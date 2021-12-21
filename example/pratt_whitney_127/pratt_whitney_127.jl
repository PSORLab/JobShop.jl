using Pkg
#=
Pkg.add(url="PATH_TO_JOBSHOP_PACKAGE HERE")
# or
=#
Pkg.rm("JobShop")
Pkg.develop(path="C:\\Users\\wilhe\\Dropbox\\My PC (DESKTOP-P6322LG)\\Desktop\\New Job Shop\\Uploaded Jobshop\\JobShop.jl")

using JobShop

jsprob = load_from_csv(@__DIR__)
jsprob.parameter.start_estimate = 1400.0
for i in 1:5
    jsprob.I[i] = Int[i]
end
for (k,i) in enumerate(6:10:106)
    jsprob.I[k + 5] = Int[j for j=i:(i+9)]
end
for (k,i) in enumerate(117:127)
    jsprob.I[k + 16] = Int[i]
end


sequential_solve!(jsprob)