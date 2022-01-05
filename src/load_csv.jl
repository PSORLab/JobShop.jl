# Load data from CSV tables into DataFrames
read_join(p,s) = CSV.read(joinpath(p,s), DataFrame, header=1)

"""
$TYPEDSIGNATURES

Load a problem formatted by five .csv files in the provided path.
"""
function load_from_csv(path)
    jsprob = JobShopProblem()
    machine              = read_join(path, "machine_capacity.csv")
    part_due             = read_join(path, "part_due.csv")
    part_operation_num   = read_join(path, "part_operation_num.csv")
    machine_part_op      = read_join(path, "machine_part_op.csv")
    part_operation_time  = read_join(path, "part_operation_time.csv")
    scrap_probability    = read_join(path, "part_operation_scrap.csv")
    rework_probability   = read_join(path, "part_operation_rework.csv")
    part_group           = read_join(path, "part_group.csv")
    
    jsprob.M = machine.capacity
    jsprob.Mi = 1:length(machine.capacity)
    jsprob.d = part_due.due
    for (kj,j) in enumerate(part_operation_num.num)
        jsprob.J[kj] = Int[i for i in 1:j]
    end

    for r in eachrow(part_group)
        if !haskey(jsprob.I, r.group)
            jsprob.I[r.group] = Int[]
        end
        push!(jsprob.I[r.group], r.part)
    end
    
    for r in eachrow(machine_part_op)
        if !haskey(jsprob.O, r.machine)
            jsprob.O[r.machine] = Tuple{Int,Int}[]
        end
        push!(jsprob.O[r.machine], (r.part,r.op))
        if !haskey(jsprob.U, (r.part,r.op))
            jsprob.U[r.part,r.op] = Int[]
        end
        push!(jsprob.U[r.part,r.op], r.machine)
    end

    for r in eachrow(part_operation_time)
        jsprob.p[r.part, r.op] = r.time
    end
    for r in eachrow(scrap_probability)
        jsprob.ps[r.part, r.op] = r.scrap
    end
    for r in eachrow(rework_probability)
        jsprob.pr[r.part, r.op] = r.rework
    end
    for m in jsprob.Mi
        append!(jsprob.Om, jsprob.O[m])
    end
    unique!(jsprob.Om)
    # TODO: use csv to specify weights

    for I ∈ values(jsprob.I), i ∈ I 
        jm1 = jsprob.J[i][1:(end-1)]
        for j = jsprob.J[i]
            jsprob.t[i,j] = -prod(k -> (1 - jsprob.ps[i,k]), jsprob.J[i])
            if !isempty(jm1)
                jsprob.t[i,j] += prod(k -> (1 - jsprob.ps[i,k]), jm1) 
            end
        end
        jsprob.ta[i] = prod(j -> (1 - jsprob.ps[i,j]),              jsprob.J[i])
        jsprob.tb[i] = sum(j -> jsprob.t[i,j]*(1 - jsprob.pr[i,j]), jsprob.J[i])
        jsprob.tc[i] = sum(j -> jsprob.t[i,j]*jsprob.pr[i,j],       jsprob.J[i])
    end

    for I ∈ values(jsprob.I), i ∈ I, j ∈ jsprob.J[i]
        jm1 = jsprob.J[i][1:(end-1)]
        if (i,j) in jsprob.Om
            pm1 = !isempty(jm1) ? prod(k -> 1 - jsprob.ps[i,k], jm1) : 0.0
            jsprob.y[i,j] = pm1
            jsprob.τ[i,j] = pm1 - prod(k -> 1 - jsprob.ps[i,k], jsprob.J[i])*jsprob.pr[i,j]
        end
    end
    jsprob.s = zeros(length(jsprob.Mi),length(jsprob.Tp))

    return jsprob
end
