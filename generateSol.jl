using JuMP
using CPLEX
using Random

include("neighborhood.jl")

function genSolWithSolver(n, m, K, B, L, l, lh, w_v, W_v, W)
    getSol = Model(CPLEX.Optimizer)

    @variable(getSol, x[e in 1:m], Int)
    @variable(getSol, z[k in 1:K, e in 1:m], Bin)
    @variable(getSol, y[k in 1:K, v in 1:n], Bin)
    @variable(getSol, alpha >= 0)
    @variable(getSol, beta[e in 1:m] >= 0)
    # from the slave problem
    @variable(getSol, zeta[k in 1:K] >= 0)
    @variable(getSol, gamma[k in 1:K, v in 1:n] >= 0)

    @constraint(getSol, [e in 1:m; node1(e,n)!=node2(e,n)], alpha + beta[e] >= x[e]*(lh[node1(e,n)]+lh[node2(e,n)]))
    @constraint(getSol, [v in 1:n], sum(y[k, v] for k in 1:K) == 1)
    @constraint(getSol, [e in 1:m, k in 1:K], y[k, node2(e,n)] + y[k, node1(e,n)] - 1 <= z[k, e])
    @constraint(getSol, [e in 1:m, k in 1:K], (y[k, node2(e,n)] + y[k, node1(e,n)]) / 2 >= z[k, e])
    @constraint(getSol, [e in 1:m], x[e] == sum(z[k, e] for k in 1:K))
    # from the slave problem
    @constraint(getSol, [k in 1:K, v in 1:n], zeta[k] + gamma[k, v] >= w_v[v]*y[k, v])
    @constraint(getSol, [k in 1:K], W*zeta[k] + sum(W_v[v]*gamma[k, v] for v in 1:n) <= B - sum(w_v[v]*y[k, v] for v in 1:n))

    @objective(getSol, Min, 0)

    # Désactive les sorties de CPLEX
    # set_optimizer_attribute(getSol, "CPX_PARAM_SCRIND", 0)
    optimize!(getSol)

    feasiblefound = primal_status(getSol) == MOI.FEASIBLE_POINT
    if feasiblefound
        clusters = JuMP.value.(y)
    end

    sol = zeros(Int64, n)
    for v in 1:n
        for k in 1:K
            if clusters[k,v]==true
                sol[v] = k
            end
        end
    end

    return sol
end

function xFromPartition(sol1D, n::Int64, m::Int64)::Vector{Bool}
    # returns a vector of bools with true if the edge is in a cluster
    return Vector([(sol1D[node1(e, n)] == sol1D[node2(e, n)]) for e in 1:m])
end

function sol2Dto1D(sol2D, K::Int64, n::Int64)
    # Input: sol is a 2D vector : sol[k][i] tells if node i is in cluster k
    # Output: vectorSol is a 1D vector : vectorSol[i] tells in which cluster node i is
    sol1D = zeros(Int64, n)
    for k in 1:K
        for i in 1:n
            if sol2D[k][i]
                sol1D[i] = k
            end
        end
    end
    return sol1D
end

function sol1Dto2D(sol1D::Array{Int64}, K::Int64, n::Int64)
    # Input: sol is a 1D vector : sol[i]=k if node i is in cluster k
    # Output: sol2D is a 2D vector : sol2D[k][i] true if node i in cluster k
    sol2D = [zeros(Bool, n) for k in 1:K]
    for v in 1:n
        sol2D[sol1D[v]][v] = true
    end
    return sol2D
end

function couplesWithSameW2(sol1D::Array{Int64}, n::Int64, w_v, w_v2)
    # returns an array of couples of vertices with compatible w_v(2) values
    couples = Vector()
    for node1 in 1:n
        cluster1 = sol1D[node1]
        for node2 in node1:n
            cluster2 = sol1D[node2]
            if cluster2!=cluster1
                if w_v[node1] < w_v2[cluster2][node2] && w_v[node2] < w_v2[cluster1][node1]
                    push!(couples, [node1, node2])
                end
            end
        end
    end
    return couples
end

function heuristic(inputFile::String, timeLimit::Int64)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    println(" ")
    println("LOOKING FOR FIRST SOLUTION")
    currentSol1D = genSolWithSolver(n, m, K, B, L, l, lh, w_v, W_v, W)
    currentSol2D = sol1Dto2D(currentSol1D, K, n)
    currentValue = partitionValue(currentSol1D, n, m, l, lh, L)
    println("found with value: ", currentValue)

    start = time()
    compTime = 0
    iter = 1

    while compTime < timeLimit
        println(" ")
        println("ITERATION ", iter)
        # trying to move around vertices
        println("*** SIMPLE MOVE ***")
        nodeIndices = [i for i in 1:n]
        shuffle!(nodeIndices)
        for nodeIdx in nodeIndices
            moved, currentValue = simpleMove(nodeIdx, currentSol1D, currentSol2D, K, B, n, m, w_v, W_v, W, l, lh, L, currentValue)
            if moved
                println("node ", nodeIdx, " moved")
            end
        end
        # currentValue = partitionValue(currentSol1D, n, m, l, lh, L)
        println("new value : ", currentValue)

        # getting the value, delta_2 and w_v2 for each cluster
        println("GETTING USEFUL DATA ABOUT SOLUTION")
        delta2 = Vector()
        w_v2 = Vector()
        clusterValues = Vector()
        for k in 1:K
            clusVal, delta2Clus = clusterValue(k, currentSol2D, n, m, w_v, W_v, W)
            push!(delta2, delta2Clus)
            append!(clusterValues, clusVal)
            w_v2k = w_v .* (delta2Clus .+ 1)
            push!(w_v2, w_v2k)
        end

        println("*** SWAPPING COUPLES ***")
        couples = couplesWithSameW2(currentSol1D, n, w_v, w_v2)
        currentValue = switchTwoNodes(currentSol1D, currentSol2D, currentValue, couples, B, n, m, w_v, W_v, W, l, lh, L)
        println("value after search : ", currentValue)

        println("**** SWAPPING NODES WITH SAME W_V ****")
        currentValue =switchNodes(currentSol1D, currentSol2D, currentValue, B, n, m, w_v, W_v, W, l, lh, L)
        currentValue = partitionValue(currentSol1D, n, m, l, lh, L)
        println("value after similar w_v swaps : ", currentValue)

        compTime = time() - start
        iter+=1
    end

    println("FINAL VALUE : ", currentValue)

    return currentSol1D, currentValue
end

# sol, val = heuristic("data/318_lin_6.tsp", 60)
