using JuMP
using CPLEX

function slave_objective(l, L::Int, l_hat, x_star, E, m::Int)
    slave_obj = Model(CPLEX.Optimizer)

    @variable(slave_obj, delta_1[e in 1:m] >= 0)

    @constraint(slave_obj, sum(delta_1[e] for e in 1:m) <= L)
    @constraint(slave_obj, [e in 1:m], delta_1[e] <= 3)

    @objective(slave_obj, Min, sum(x_star[e]*(l[e]+delta_1[e]*(l_hat[E[e][1]]+l_hat[E[e][2]])) for e in 1:m))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(slave_obj, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(slave_obj, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(slave_obj, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    # set_optimizer_attribute(slave_obj, "CPX_PARAM_SCRIND", 0)

    # start = time()
    # optimize!(slave_obj)
    # computation_time = time() - start

    feasiblefound = primal_status(slave_obj) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(slave_obj)
    end

    return obj# , computation_time
    # SHALL RETURN DELTA1_STAR
end

function obj_sub_pb(l, L::Int, l_hat, x_star, E, m::Int)
    

function slave_constraint(k::Int, w, W_global::Int, W, y_star, V, n::Int)
    slave_cons = Model(CPLEX.Optimizer)

    @variable(slave_cons, delta_2[v in 1:n] >= 0)

    @constraint(slave_cons, sum(delta_2[v] for v in 1:v) <= W_global)
    @constraint(slave_cons, [v in 1:n], delta_2[v] <= W[v])

    @objective(slave_cons, Max, sum(y_star[k][v]*w[v]*(1+delta_2[v]) for v in 1:v))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(slave_cons, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(slave_cons, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(slave_cons, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    # set_optimizer_attribute(slave_cons, "CPX_PARAM_SCRIND", 0)

    # start = time()
    # optimize!(slave_cons)
    # computation_time = time() - start

    feasiblefound = primal_status(slave_cons) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(slave_cons)
    end

    return obj# , computation_time
    # SHALL RETURN DELTA2_STAR
end
