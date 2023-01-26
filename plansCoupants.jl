include("coupes.jl")


function master_pb(n::Int, m::Int, K::Int, w_v, W::Int, W_v, l, L::Int, lh, B::Int, U_1, U_2)
    master = Model(CPLEX.Optimizer)

    @variable(master, x[e in 1:m], Int)
    @variable(master, z[k in 1:K, e in 1:m], Bin)
    @variable(master, y[k in 1:K, v in 1:n], Bin)
    @variable(master, t >= 0)

    @constraint(master, [v in 1:n], sum(y[k, v] for k in 1:K) == 1)
    @constraint(master, [e in 1:m, k in 1:K], y[k, node1(e, n)] + y[k, node2(e, n)] - 1 <= z[k, e])
    @constraint(master, [e in 1:m, k in 1:K], ((y[k, node1(e, n)] + y[k, node2(e, n)])) / 2 >= z[k, e])
    @constraint(master, [e in 1:m], x[e] == sum(z[k, e] for k in 1:K))
    # slave objective
    @constraint(master, [i in 1:length(U_1)], t >= sum(U_1[i][e]*x[e] for e in 1:m))
    # slave constraint
    @constraint(master, [i in 1:length(U_2), k in 1:K], sum(U_2[i][v]*y[k, v] for v in 1:n) <= B)

    @objective(master, Min, t)

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(master, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0)

    # start = time()
    optimize!(master)
    # computation_time = time() - start

    feasiblefound = primal_status(master) == MOI.FEASIBLE_POINT
    if feasiblefound
        t_star = JuMP.value.(t)
        x_star = JuMP.value.(x)
        y_star = JuMP.value.(y)
        obj = JuMP.objective_value(master)
    end

    return t_star, x_star, y_star, obj #, computation_time
end

function solveByCuts(inputFile::String, timeLimit::Int64)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W_v, K, B
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    U_1 = Vector{Vector{Float64}}()
    push!(U_1, [l[e]/2 for e in 1:m])
    U_2 = Vector{Vector{Float64}}()
    push!(U_2, [w_v[v] for v in 1:n])

    start = time()

    # first solution
    continue_resolution = false
    println("FIRST MASTER")
    t, x, y, lower_bound = master_pb(n, m, K, w_v, W, W_v, l, L, lh, B, U_1, U_2)
    println("FIRST SLAVE")
    delta1, obj_val = slave_objective(l, L, lh, x, m)
    if obj_val > t
        new_constraint = [1/2*l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)]) for e in 1:m]
        push!(U_1, new_constraint)
        continue_resolution = true
    end
    for k in 1:K
        delta2, cons_val = slave_constraint(k, w_v, W, W_v, y, n)
        if cons_val > B
            new_constraint = [w_v[v]*(1+delta2[v]) for v in 1:n]
            push!(U_2, new_constraint)
            continue_resolution = true
        end
    end

    # iterations of the resolution
    while continue_resolution && time()-start<timeLimit
        println("TIME: ", time()-start)
        continue_resolution = false
        t, x, y, lower_bound = master_pb(n, m, K, w_v, W, W_v, l, L, lh, B, U_1, U_2)
        delta1, obj_val = slave_objective(l, L, lh, x, m)
        if obj_val > t
            new_constraint = [1/2*l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)]) for e in 1:m]
            push!(U_1, new_constraint)
            continue_resolution = true
        end
        for k in 1:K
            delta2, cons_val = slave_constraint(k, w_v, W, W_v, y, n)
            if cons_val > B
                new_constraint = [w_v[v]*(1+delta2[v]) for v in 1:n]
                push!(U_2, new_constraint)
                continue_resolution = true
            end
        end
    end


    computation_time = time() - start

    if continue_resolution
        println("time limit reached")
        status = "time limit"
    else
        if t==lower_bound
            println("found solution hitting lower bound")
            status = "lower bound"
        else
            println("solution found because constraints not violated")
            status = "valid solution"
        end
    end

    res = [[] for k in 1:K]
    for k in 1:K
        for v in 1:n
            if y[k,v]==true
                append!(res[k], v)
            end
        end
    end

    return res, t, lower_bound, computation_time, status
end


# inputFile = "data/10_ulysses_3.tsp"
# res, t, lower_bound, computation_time = solveByCuts(inputFile)
# println("clusters : ", res)
# println("value : ", t)
# println("computation time : ", computation_time)
