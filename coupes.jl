using JuMP
using CPLEX

include("common.jl")

function slave_objective(l, L::Int, lh, x_star, m::Int)
    slave_obj = Model(CPLEX.Optimizer)

    @variable(slave_obj, delta_1[e in 1:m] >= 0)

    @constraint(slave_obj, [i in 1:n], delta_1[getEdge(n, i, i)] == 0)
    @constraint(slave_obj, sum(delta_1[e] for e in 1:m) <= L)
    @constraint(slave_obj, [e in 1:m], delta_1[e] <= 3)

    @objective(slave_obj, Max, sum(x_star[e]*(1/2*l[e]+delta_1[e]*(lh[node1(e, n)]+lh[node2(e, n)])) for e in 1:m))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(slave_obj, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(slave_obj, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(slave_obj, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(slave_obj, "CPX_PARAM_SCRIND", 0)

    # start = time()
    optimize!(slave_obj)
    # computation_time = time() - start

    feasiblefound = primal_status(slave_obj) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta1_star = JuMP.value.(delta_1)
        obj = JuMP.objective_value(slave_obj)
    end

    return delta1_star, obj #, computation_time
    # SHALL RETURN DELTA1_STAR
end

function slave_constraint(k::Int, w_v, W::Int, W_v, y_star, n::Int)
    slave_cons = Model(CPLEX.Optimizer)

    @variable(slave_cons, delta_2[v in 1:n] >= 0)

    @constraint(slave_cons, sum(delta_2[v] for v in 1:n) <= W)
    @constraint(slave_cons, [v in 1:n], delta_2[v] <= W_v[v])

    @objective(slave_cons, Max, sum(y_star[k,v]*w_v[v]*(1+delta_2[v]) for v in 1:n))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(slave_cons, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(slave_cons, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(slave_cons, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(slave_cons, "CPX_PARAM_SCRIND", 0)

    # start = time()
    optimize!(slave_cons)
    # computation_time = time() - start

    feasiblefound = primal_status(slave_cons) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta2_star = JuMP.value.(delta_2)
        obj = JuMP.objective_value(slave_cons)
    end

    return delta2_star, obj #, computation_time
    # SHALL RETURN DELTA2_STAR
end

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

function solveByCuts(inputFile::String)
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

    # iterations of the resolution
    while continue_resolution
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

    if t==lower_bound
        println("found solution hitting lower bound")
    else
        println("solution found because constraints not violated")
    end


    res = [[] for k in 1:K]
    for k in 1:K
        for v in 1:n
            if y[k,v]==true
                append!(res[k], v)
            end
        end
    end

    return res, t, lower_bound, computation_time
end


function master_pb_via_CPLEX(inputFile::String)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W_v, K, B
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    U_1 = [l[e] for e in 1:m]
    U_2 = [w_v[v] for v in 1:n]

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
    @constraint(master, t >= sum(U_1[e]*x[e] for e in 1:m))
    # slave constraint
    @constraint(master, [k in 1:K], sum(U_2[v]*y[k, v] for v in 1:n) <= B)

    @objective(master, Min, t)

    # Désactive le presolve (simplification automatique du modèle)
    set_optimizer_attribute(master, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    set_optimizer_attribute(master, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0)

    callback_called = false
    function callback_function(cb_data::CPLEX.CallbackContext)
        callback_called = true
        # println("USING THE CALLBACK FUNCTION YEAY")
        # add in arguments context_id::Clong ?
        # CPLEX.load_callback_variable_primal(cb_data, context_id)

        t_star = callback_value.(Ref(cb_data), t)
        x_star = callback_value.(Ref(cb_data), x)
        y_star = callback_value.(Ref(cb_data), y)

        delta1, obj_val = slave_objective(l, L, lh, x_star, m)
        if obj_val > t_star
            println(" !!!!!!!!!!!!!!! adding objective constraint")
            new_constraint = [l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)]) for e in 1:m]
            cut = @build_constraint(sum(x[e]*new_constraint[e] for e in 1:m) <= t)
            # print("NEW CONSTRAINT : ", cut)
            MOI.submit(master, MOI.UserCut(cb_data), cut)
        end

        for k in 1:K
            delta2, cons_val = slave_constraint(k, w_v, W, W_v, y_star, n)
            if cons_val > B
                println(" !!!!!!!!!!!!!!! adding constraint in k = ", k)
                cut = @build_constraint(sum(y[v]*w_v[v]*(1+delta2[v]) for v in 1:n) <= B)
                # print("NEW CONSTRAINT : ", cut)
                MOI.submit(master, MOI.UserCut(cb_data), cut)
            end
        end

    end

    MOI.set(master, MOI.UserCutCallback(), callback_function)

    start = time()
    optimize!(master)
    computation_time = time() - start

    @show callback_called
    println()
    # print(master)

    feasible_found = primal_status(master) == MOI.FEASIBLE_POINT
    if !feasible_found
        println("No feasible point found !")
        return
    end

    t_star = JuMP.value.(t)
    x_star = JuMP.value.(x)
    y_star = JuMP.value.(y)
    obj = JuMP.objective_value(master)

    res = [[] for k in 1:K]
    for k in 1:K
        for v in 1:n
            if y_star[k,v]==true
                append!(res[k], v)
            end
        end
    end

    return res, t_star, obj, computation_time
end

inputFile = "data/10_ulysses_3.tsp"
clusters, t_star, value, computation_time = master_pb_via_CPLEX(inputFile)
println("clusters : ", clusters)
println("value : ", value)
println("computation time : ", computation_time)
