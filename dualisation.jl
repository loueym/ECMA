using JuMP
using CPLEX

function master_pb(inputFile::String)
    include(inputFile)          # contains G=(V,E), l[e], l_hat[v], L, w[v], W[v], W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    n = length(V)
    m = length(E)

    master = Model(CPLEX.Optimizer)

    @variable(master, x[e in 1:m], Int)
    @variable(master, z[k in 1:K, e in 1:m], Bin)
    @variable(master, y[k in 1:K, v in 1:n], Bin)
    @variable(master, alpha >= 0)
    @variable(master, beta[e in 1:m] >= 0)
    @variable(master, t[k in 1:K]) # for the slave problem


    @constraint(master, [e in 1:m], alpha + beta[e] >= x[e]*(l_hat[E[e][1]]+l_hat[E[e][2]]))
    @constraint(master, [v in 1:n], sum(y[k][v] for k in 1:K) == 1)
    @constraint(master, [e in 1:m, k in 1:K], sum(y[k][v]*A[v][e] for v in 1:n) - 1 <= z[k][e])
    @constraint(master, [e in 1:m, k in 1:K], sum(y[k][v]*A[v][e] for v in 1:n) / 2 >= z[k][e])
    @constraint(master, [e in 1:m], x[e] == sum(z[k][e] for k in 1:K))
    # slave problem
    @constraint(master, [k in 1:K], t[k] == slave_pb(k, n, w, y, W))
    @constraint(master, [k in 1:K], t[k] <= B - sum(w[v]*y[k][v] for v in 1:n))

    @objective(master, Min, sum(x[e]*l[e] for e in 1:m) + L*alpha + 3*sum(beta[e] for e in 1:m))

    # Désactive le presolve (simplification automatique du modèle)
    set_optimizer_attribute(master, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    set_optimizer_attribute(master, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    set_optimizer_attribute(master, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0)

    start = time()
    optimize!(master)
    computation_time = time() - start

    feasiblefound = primal_status(master) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(master)
    end

    return obj, computation_time
end

function slave_pb(k::Int, n::Int, w, y, W_global::Int, W)
    slave = Model(CPLEX.Optimizer)

    @variable(slave, zeta >= 0)
    @variable(slave, gamma[v in 1:n] >= 0)

    @constraint(slave, [v in 1:n], zeta + gamma[v] >= w[v]*y[k][v])

    @objective(slave, Min, W_global*zeta + sum(W[v]*gamma[v] for v in 1:n))

    # Désactive le presolve (simplification automatique du modèle)
    set_optimizer_attribute(slave, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    set_optimizer_attribute(slave, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    set_optimizer_attribute(slave, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(slave, "CPX_PARAM_SCRIND", 0)

    # start = time()
    # optimize!(slave)
    # computation_time = time() - start

    feasiblefound = primal_status(slave) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(slave)
    end

    return obj# , computation_time
end
