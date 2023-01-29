using JuMP
using CPLEX

include("common.jl")
include("neighborhood.jl")

function xFromPartition(vectorPartition, n::Int64, m::Int64)::Vector{Bool}
    # returns a vector of bools with true if the edge is in a cluster
    return Vector([(vectorPartition[node1(e, n)] == vectorPartition[node2(e, n)]) for e in 1:m])
end

function partitionValue(vectorPartition, n::Int64, m::Int64, l, lh, L::Int64)::Float64
    # returns the value of the given partition
    x = xFromPartition(vectorPartition, n, m)

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

function clusterValue(k::Int64, partition, n::Int64, m::Int64, w_v, W_v, W::Int64)
    # rerturns the value of cluster k in the given partition
    clusVal = Model(CPLEX.Optimizer)
    @variable(clusVal, delta2[v in 1:n] >= 0)

    @constraint(clusVal, sum(delta2[v] for v in 1:n) <= W)
    @constraint(clusVal, [v in 1:n], delta2[v] <= W_v[v])

    @objective(clusVal, Max, sum(w_v[v]*partition[k][v]*(1+delta2[v]) for v in 1:n))
    optimize!(clusVal)

    feasiblefound = primal_status(clusVal) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta2Star = JuMP.value.(delta2)
        clusterVal = JuMP.objective_value(clusVal)
    end

    return clusterVal, delta2Star
end

function genRandomSol(n::Int64, m::Int64, K::Int64, B::Int64, L::Int64, l, lh, w_v, W_v, W::Int64)
    # returns a randomly generated admissible solution
    sol = [[false for i in 1:n] for k in 1:K]
    nodes = [i for i in 1:n]
    for k in 1:K
        idx = rand(1:length(nodes))
        sol[k][nodes[idx]] = true
        deleteat!(nodes, idx)
    end
    iter = 0
    delta2Star = nothing
    while length(nodes)>0 && iter < n+2*K
        # picking a cluster to extend
        clusterIdx = rand(1:K)
        # taking a random node and trying to add it to the cluster
        idx = rand(1:length(nodes))
        sol[clusterIdx][nodes[idx]] = true
        # testing the value of the cluster
        clusterVal, delta2Star = clusterValue(clusterIdx, sol, n, m, w_v, W_v, W)
        if clusterVal <= B
            deleteat!(nodes, idx)
        else
            sol[clusterIdx][nodes[idx]] = false
        end
        iter += 1
    end
    return sol, delta2Star
end

function couplesWithSameW2(w_v2)
    # returns an array of couples of vertices with same w_2 value
    n = length(w_v2)
    couples = []
    for d in Set(w_v2)
        indices = findall(x->x==d, w_v2)
        m = length(indices)
        indexInVector = 1
        for i in indices
            for j in indices[indexInVector+1:m]
                push!(couples, (i, j))
            end
            indexInVector += 1
        end
    end
    return couples
end

function sol2Dto1D(sol)
    # Input: sol is a 2D vector : sol[k][i] tells if node i is in cluster k
    # Output: vectorSol is a 1D vector : vectorSol[i] tells in which cluster node i is
    K = length(sol)
    n = length(sol[1])
    vectorSol = [0 for i in 1:n]
    for k in 1:K
        for i in 1:n
            if sol[k][i]
                vectorSol[i] = k
            end
        end
    end
    return vectorSol
end


function heuristic(inputFile::String, timeLimit::Int64)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    # does not always generate a full solution
    sol, delta2 = genRandomSol(n, m, K, B, L, l, lh, w_v, W_v, W)
    vectorSol = sol2Dto1D(sol)
    val = partitionValue(vectorSol, n, m, l, lh, L)
    println("First solution value: ", val)

    w_v2 = w_v .* (delta2 .+ 1)

    # generating a list of tuples of nodes with same w_v^2
    couples = couplesWithSameW2(w_v2)

    start = time()
    computationTime = time() - start
    while computationTime < timeLimit
        # we try to switch two nodes (which are not currently in the same cluster)
        vectorSol2 = switchTwoNodes(vectorSol, couples)
        val2 = partitionValue(vectorSol2, n, m, l, lh, L)
        if val2 < val
            val = val2
            vectorSol = vectorSol2
            println("After ", computationTime, " seconds, found a better solution with value ", val)
        end

        computationTime = time() - start
    end

    println("generated sol : ", vectorSol)
    return vectorSol
end
