using JuMP
using CPLEX

include("common.jl")

function master_pb(inputFile::String)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n

    l = generate_l(coordinates, n)

    master = Model(CPLEX.Optimizer)

    @variable(master, x[e in 1:m], Int)
    @variable(master, z[k in 1:K, e in 1:m], Bin)
    @variable(master, y[k in 1:K, v in 1:n], Bin)
    @variable(master, alpha >= 0)
    @variable(master, beta[e in 1:m] >= 0)
    # @variable(master, t[k in 1:K]) # for the slave problem
    # from the slave problem
    @variable(master, zeta[k in 1:K] >= 0)
    @variable(master, gamma[k in 1:K, v in 1:n] >= 0)


    @constraint(master, [e in 1:m; node1(e,n)!=node2(e,n)], alpha + beta[e] >= x[e]*(lh[node1(e,n)]+lh[node2(e,n)]))
    @constraint(master, [v in 1:n], sum(y[k, v] for k in 1:K) == 1)
    @constraint(master, [e in 1:m, k in 1:K], y[k, node2(e,n)] + y[k, node1(e,n)] - 1 <= z[k, e])
    @constraint(master, [e in 1:m, k in 1:K], (y[k, node2(e,n)] + y[k, node1(e,n)]) / 2 >= z[k, e])
    @constraint(master, [e in 1:m], x[e] == sum(z[k, e] for k in 1:K))
    # slave problem
    #@constraint(master, [k in 1:K], t[k] == slave_pb(k, n, y, W, W_v, w_v))
    #@constraint(master, [k in 1:K], t[k] <= B - sum(w_v[v]*y[k, v] for v in 1:n))
    # from the slave problem
    @constraint(master, [k in 1:K, v in 1:n], zeta[k] + gamma[k, v] >= w_v[v]*y[k, v])
    @constraint(master, [k in 1:K], W*zeta[k] + sum(W_v[v]*gamma[k, v] for v in 1:n) <= B - sum(w_v[v]*y[k, v] for v in 1:n))

    @objective(master, Min, sum(x[e]*l[e]/2 for e in 1:m) + L*alpha + 3*sum(beta[e] for e in 1:m))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(master, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(master, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    # set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0)

    start = time()
    optimize!(master)
    computation_time = time() - start

    feasiblefound = primal_status(master) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(master)
        clusters = JuMP.value.(y)
    end

    res = [[] for k in 1:K]
    for k in 1:K
        for v in 1:n
            if clusters[k,v]==true
                append!(res[k], v)
            end
        end
    end
    xStar = [false for i in 1:m]
    for k in 1:K
        for v1 in res[k]
            for v2 in res[k]
                xStar[(v1-1)*n + v2] = true
            end
        end
    end

    return xStar, res, obj, computation_time
end

# INUTILE (intégré dans le master_pb)
function slave_pb(k::Int, n::Int, y, W::Int, W_v, w_v)
    slave = Model(CPLEX.Optimizer)

    @variable(slave, zeta >= 0)
    @variable(slave, gamma[v in 1:n] >= 0)

    @constraint(slave, [v in 1:n], zeta + gamma[v] >= w_v[v]*y[k, v])

    @objective(slave, Min, W*zeta + sum(W_v[v]*gamma[v] for v in 1:n))

    # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(slave, "CPXPARAM_Preprocessing_Presolve", 0)
    # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(slave, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # Désactive la génération de solutions entières à partir de solutions fractionnaires
    # set_optimizer_attribute(slave, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # Désactive les sorties de CPLEX (optionnel)
    # set_optimizer_attribute(slave, "CPX_PARAM_SCRIND", 0)

    # start = time()
    # optimize!(slave)
    # computation_time = time() - start

    feasiblefound = primal_status(slave) == MOI.FEASIBLE_POINT
    if feasiblefound
        obj = JuMP.objective_value(slave)
    end

    return obj# , computation_time
end

xStar, clustersTest, test, testCompTime = master_pb("data/14_burma_3.tsp")
println("attained value: ", test)
println("time needed: ", testCompTime)
println("clusters: ", clustersTest)
