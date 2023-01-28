include("coupes.jl")

function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    # Fonction permettant de déterminer si c’est l’obtention d’une
    # solution entière qui a entraîné l’appel d’un callback
    # (il n’est pas nécessaire d’en comprendre le fonctionnement)
    # context_id == CPX_CALLBACKCONTEXT_CANDIDATE si le callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    # S’il n’y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

function solveByBnC(inputFile::String, timeLimit::Int64)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W_v, K, B
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    master = Model(CPLEX.Optimizer)
    MOI.set(master, MOI.NumberOfThreads(), 1)

    @variable(master, x[e in 1:m], Int)
    @variable(master, z[k in 1:K, e in 1:m], Bin)
    @variable(master, y[k in 1:K, v in 1:n], Bin)
    @variable(master, t >= 0)

    @constraint(master, [v in 1:n], sum(y[k, v] for k in 1:K) == 1)
    @constraint(master, [e in 1:m, k in 1:K], y[k, node1(e, n)] + y[k, node2(e, n)] - 1 <= z[k, e])
    @constraint(master, [e in 1:m, k in 1:K], ((y[k, node1(e, n)] + y[k, node2(e, n)])) / 2 >= z[k, e])
    @constraint(master, [e in 1:m], x[e] == sum(z[k, e] for k in 1:K))
    # slave objective
    @constraint(master, t >= sum(1/2*l[e]*x[e] for e in 1:m))
    # slave constraint
    @constraint(master, [k in 1:K], sum(w_v[v]*y[k, v] for v in 1:n) <= B)

    @objective(master, Min, t)

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(master, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0)
    set_time_limit_sec(master, timeLimit)

    callback_called = false
    function callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
            callback_called = true
            # println("USING THE CALLBACK FUNCTION YEAY")
            # add in arguments context_id::Clong ?
            CPLEX.load_callback_variable_primal(cb_data, context_id)

            x_star = callback_value.(cb_data, x)
            t_star = callback_value.(cb_data, t)
            y_star = callback_value.(cb_data, y)

            # t_star = callback_value.(Ref(cb_data), t)
            # x_star = callback_value.(Ref(cb_data), x)
            # y_star = callback_value.(Ref(cb_data), y)

            delta1, obj_val = slave_objective(l, L, lh, x_star, m)
            if obj_val > t_star # - 1e-5
                #println(" !!!!!!!!!!!!!!! adding objective constraint because obj value is ", obj_val, " and t is ", t_star)
                #sleep(0.1)
                new_constraint = [1/2*l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)]) for e in 1:m]
                cut = @build_constraint(sum(x[e]*new_constraint[e] for e in 1:m) <= t)
                # print("NEW CONSTRAINT : ", cut)
                MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
            end

            for k in 1:K
                delta2, cons_val = slave_constraint(k, w_v, W, W_v, y_star, n)
                if cons_val > B # - 1e-5
                    #println(" !!!!!!!!!!!!!!! adding constraint in k = ", k, " because constraint_value is ", cons_val, " and b is ", B)
                    #sleep(0.1)
                    cut = @build_constraint(sum(y[k, v]*w_v[v]*(1+delta2[v]) for v in 1:n) <= B)
                    # print("NEW CONSTRAINT : ", cut)
                    MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
                end
            end

        end
    end

    MOI.set(master, CPLEX.CallbackFunction(), callback_function)

    start = time()
    optimize!(master)
    computation_time = time() - start

    @show callback_called
    # print(master)
    gap = 1 - MOI.get(master, MOI.RelativeGap())
    println("gap info ", gap)
    println()

    feasible_found = primal_status(master) == MOI.FEASIBLE_POINT
    if !feasible_found
        println("No feasible point found !")
        return [], [], 0, computation_time, 1, "time limit"
    end

    t_star = JuMP.value.(t)
    x_star = JuMP.value.(x)
    y_star = JuMP.value.(y)
    obj = JuMP.objective_value(master)

    # for c in all_constraints(master, AffExpr, MOI.LessThan{Float64})
    #     println(c)
    # end

    res = [[] for k in 1:K]
    for k in 1:K
        for v in 1:n
            if y_star[k,v]==true
                append!(res[k], v)
            end
        end
    end

    return res, t_star, obj, computation_time, gap, "solved"
end


# inputFile = "data/38_rat_6.tsp"
# clusters, t_star, value, computation_time, gap, status = solveByBnC(inputFile, 30)
# println("status: ", status)
# println("clusters : ", clusters)
# println("value : ", value)
# println("computation time : ", computation_time)
