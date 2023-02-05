### SUB-PROBLEMS FOR BnC AND CUTS METHODS

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

    optimize!(slave_obj)

    feasiblefound = primal_status(slave_obj) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta1_star = JuMP.value.(delta_1)
        obj = JuMP.objective_value(slave_obj)
    end

    return delta1_star, obj
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

    optimize!(slave_cons)

    feasiblefound = primal_status(slave_cons) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta2_star = JuMP.value.(delta_2)
        obj = JuMP.objective_value(slave_cons)
    end

    return delta2_star, obj
end
