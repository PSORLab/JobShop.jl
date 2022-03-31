module JobShop

using CSV, DataFrames, DocStringExtensions, JuMP, Printf, Requires, UnPack

export configure!, load_from_csv, sequential_solve!, Ext, JobShopProblem

include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "load_csv.jl"))
include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "display.jl"))
include(joinpath(@__DIR__, "feasibility_problem.jl"))
include(joinpath(@__DIR__, "subproblem.jl"))
include(joinpath(@__DIR__, "sequential_solve.jl"))
#include(joinpath(@__DIR__, "parallel_solve.jl"))

function __init__()
    @require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0"      include(joinpath(@__DIR__, "config", "cplex.jl"))
    @require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b"     include(joinpath(@__DIR__, "config", "gurobi.jl"))
    #@require KNITRO="67920dd8-b58e-52a8-8622-53c4cffbe346"     include(joinpath(@__DIR__, "config", "knitro.jl"))
    #@require MosekTools="1ec41992-ff65-5c91-ac43-2df89e9693a4" include(joinpath(@__DIR__, "config", "mosek.jl"))
    #@require Xpress="9e70acf3-d6c9-5be6-b5bd-4e2c73e3e054"     include(joinpath(@__DIR__, "config", "xpress.jl"))
end

end
