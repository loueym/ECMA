using JuMP
using CPLEX

include("common.jl")

function partitionValue(partition, n, m, l, lh, L)
    solVal = Model(CPLEX.Optimizer)
    @variable(solVal, delta1[e in 1:m] >= 0)

    @constraint(solVal, sum(delta1[e] for e in 1:m) <= L)
    @constraint(solVal, [e in 1:m], delta1[e] <= 3)

    @objective(solVal, Max, sum(x[e]*(1/2*l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)])) for e in 1:m))
    optimize!(solVal)

    feasiblefound = primal_status(solVal) == MOI.FEASIBLE_POINT
    if feasiblefound
        value = JuMP.objective_value(solVal)
    end

    return value
end

function clusterValue(k, partition, n, m, w_v, W_v, W)
    clusVal = Model(CPLEX.Optimizer)
    @variable(clusVal, delta2[v in 1:n] >= 0)

    @constraint(clusVal, sum(delta2[v] for v in 1:n) <= W)
    @constraint(clusVal, [v in 1:n], delta2[v] <= W_v[v])

    @objective(clusVal, Max, sum(w_v[v]*partition[k][v]*(1+delta2[v]) for v in 1:n))
    optimize!(clusVal)

    feasiblefound = primal_status(clusVal) == MOI.FEASIBLE_POINT
    if feasiblefound
        clusterVal = JuMP.objective_value(clusVal)
    end

    return clusterVal
end

function genRandomSol(n, m, K, B, L, l, lh, w_v, W_v, W)
    sol = [[false for i in 1:n] for k in 1:K]
    nodes = [i for i in 1:n]
    for k in 1:K
        idx = rand(1:length(nodes))
        sol[k][nodes[idx]] = true
        deleteat!(nodes, idx)
    end
    iter = 0
    while length(nodes)>0 && iter < n+2*K
        # picking a cluster to extend
        clusterIdx = rand(1:K)
        # a node to add
        idx = rand(1:length(nodes))
        sol[clusterIdx][nodes[idx]] = true
        # testing the value of the cluster
        clusterVal = clusterValue(clusterIdx, sol, n, m, w_v, W_v, W)
        if clusterVal <= B
            deleteat!(nodes, idx)
        else
            sol[clusterIdx][nodes[idx]] = false
        end
        iter += 1
    end
    return sol
end



function heuristic(inputFile::String)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    # does not always generate a full solution
    sol = genRandomSol(n, m, K, B, L, l, lh, w_v, W_v, W)

    println("generated sol : ", sol)
end
