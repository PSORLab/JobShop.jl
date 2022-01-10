# Load data from CSV tables into DataFrames
read_join(p::String,s::String) = CSV.read(joinpath(p,s), DataFrame, header=1)

"""
$TYPEDSIGNATURES

Load a problem defined by .csv files in the provided path.
"""
function load_from_csv(path::String)
    jsp = JobShopProblem()
    
    part_due           = read_join(path, "part_due.csv")
    machine_part_op    = read_join(path, "machine_part_op.csv")
    machine            = read_join(path, "machine_capacity.csv")
    part_operation_num = read_join(path, "part_operation_num.csv")
    PartOpTs           = read_join(path, "part_operation_time.csv")
    part_group         = read_join(path, "part_group.csv")   

    MaPartOps = machine_part_op
    jsp.PartDue    = part_due.due
    jsp.MachineCap = machine.capacity
    jsp.MachineType = 1:length(machine.capacity)
    for r in eachrow(part_group)
        if !haskey(jsp.Ii, r.group)
            jsp.Ii[r.group] = Int[]
        end
        push!(jsp.Ii[r.group], r.part)
    end
    for i in keys(jsp.Ii)
        append!(jsp.I, jsp.Ii[i])
    end
    unique!(jsp.I)
    
    for (kj,j) in enumerate(part_operation_num.num)
        jsp.J[kj] = j
        jsp.Jop[kj] = 1:j
    end

    nbPart = length(jsp.I)
    nbOperation = length(jsp.J)
    nbMachine = length(jsp.MachineType)

    Pro = zeros(nbPart, nbOperation)
    APro = zeros(nbPart, nbOperation)
    BPro = zeros(nbPart, nbOperation)
    for r in eachrow(PartOpTs)
        Pro[r.part,r.op] = r.time
    end
    for i in jsp.I, j in jsp.J[i]
        if j <= jsp.J[i] 
            for l in jsp.J[i]            
                if l <= j 
                    APro[i,j] += Pro[i,l]
                end
            end
        end
    end
    for i in jsp.I, j in jsp.J[i]
        if  j <= jsp.J[i]            
            BPro[i,j] = APro[i,jsp.J[i]] - APro[i,j]
        end
    end
    PartOpsMa = Tuple{Int,Int,Int,Int}[]                                                 
    for m1 in eachrow(MaPartOps), m2 in eachrow(MaPartOps)
        if (m1.part == m2.part) && (m1.op == m2.op) && (m1.machine <= m2.machine-1) 
            push!(PartOpsMa, (m1.part, m1.op, m1.machine, m2.machine))
        end
    end
    for P in eachrow(PartOpTs)
        push!(jsp.IJT, PartOpT(P.part, P.op, P.time))
    end
    for P in eachrow(MaPartOps)
        push!(jsp.MIJ, MaPartOp(P.machine, P.part, P.op))
    end

    for i in jsp.I, j in jsp.Jop[i], t in jsp.T
        jsp.sbI1[i,j,t] = 0.0 # 0.0
        for j1 in jsp.Jop[i], r in jsp.R
            jsp.sbI2[i,j,j1,r,t] = 0.0 #1.0 # 0.0
        end
    end
    
    jsp.sslackk = zeros(Float64,nbMachine,length(jsp.T))
    jsp.sv_p = zeros(Float64,nbMachine,length(jsp.T))
    jsp.mult = zeros(Float64,nbMachine,length(jsp.T))

    jsp.sTard1 = zeros(nbPart)
    jsp.sTard2 = zeros(nbPart, nbOperation, length(jsp.R)) 
    jsp.sbTime1 = zeros(nbPart,nbOperation)
    return jsp
end

#=
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
=#