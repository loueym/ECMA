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
    set_optimizer_attribute(getSol, "CPX_PARAM_SCRIND", 0)
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

function xFromPartition(vectorPartition, n::Int64, m::Int64)::Vector{Bool}
    # returns a vector of bools with true if the edge is in a cluster
    return Vector([(vectorPartition[node1(e, n)] == vectorPartition[node2(e, n)]) for e in 1:m])
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

### TO GENERATE A RANDOM SOLUTION (DOES NOT ALWAYS SUCCEED)
function genRandomSol(n::Int64, m::Int64, K::Int64, B::Int64, L::Int64, l, lh, w_v, W_v, W::Int64)
    # returns a randomly generated admissible solution
    sol = [[false for i in 1:n] for k in 1:K]
    nodes = [i for i in 1:n]
    for k in 1:K
        idx = rand(1:length(nodes))
        sol[k][nodes[idx]] = true
        deleteat!(nodes, idx)
    end
    println("after first round: ", sol)

    iter = 0
    delta2Star = nothing
    while length(nodes)>0 && iter < n+2*K
        # picking a cluster to extend
        clusterIdx = rand(1:K)
        # taking a random node and trying to add it to the cluster
        idx = rand(1:length(nodes))
        println("trying node : ", nodes[idx])
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

    if length(nodes)==0
        println("STOP bcs : all nodes dispatched")
        solFound = true
    else
        println("STOP bcs : too many iterations")
        solFound = false
    end
    return solFound, sol, delta2Star
end

function affectNodes(sol, nodes::Array{Int64}, iterLimit::Float64, k::Int64, direction::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64, B::Int64)
    iter = 0
    while isempty(nodes)==false && iter < iterLimit
        # taking the first remaining node, to add it to the cluster
        sol[k][nodes[1]] = true
        # testing the value of the cluster
        clusterVal, delta2Star = clusterValue(k, sol, n, m, w_v, W_v, W)
        if clusterVal <= B
            deleteat!(nodes, 1)
        else
            sol[k][nodes[1]] = false
            buffer = nodes[1]
            deleteat!(nodes, 1)
            push!(nodes, buffer)
        end
        # updating cluster to extend
        k += direction
        if k==1 || k==K
            direction *= -1
        end
        iter += 1
    end
end

function affectFirstNodes(n::Int64, m::Int64, K::Int64, B::Int64, L::Int64, l, lh, w_v, W_v, W::Int64)
    sol = [[false for i in 1:n] for k in 1:K]
    nodes = [i for i in 1:n]
    sort!(nodes, by=x->w_v[x], rev=true) # nodes ordered by descending static weight
    println("nodes: ", nodes)
    # a heavy node per cluster
    k = 1
    direction = 1
    iter = 0
    while iter < K
        sol[k][nodes[1]] = true
        deleteat!(nodes, 1)
        k += direction
        if k==1 || k==K
            direction *= -1
        end
        iter += 1
    end

    # affect other nodes if possible
    reverse!(nodes) # nodes ordered by ascending static weight
    affectNodes(sol, nodes, n*11/10, k, direction, n, m, w_v, W_v, W, B)

    # check if all nodes are affected
    if isempty(nodes)
        println("ALL NODES DISPATCHED")
    else
        println("TOO MANY ITERATIONS")
    end
    return isempty(nodes), sol, nodes
end

function affectNode(nodeIdx::Int64, sol, K::Int64, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64)
    k = 1 # cluster to put the node in
    goOn = true
    while k<=K && goOn
        # pick the node to move
        nodesInCluster = Vector{Int64}()
        for idx in 1:n
            if sol[k][idx]
                append!(nodesInCluster, idx)
            end
        end
        nodeToMove = argmax(x->w_v[x], nodesInCluster)
        clusterToTest = mod(k+1, K) # cluster to move nodeToMove to
        if clusterToTest==0
            clusterToTest = K
        end
        # test where to move it to
        possibleSol = copy(sol)
        tryCluster = true
        while clusterToTest!=k && tryCluster
            possibleSol[clusterToTest][nodeToMove] = true
            movedValue, delta2 = clusterValue(clusterToTest, possibleSol, n, m, w_v, W_v, W)
            if movedValue <= B
                possibleSol[k][nodeToMove] = false
                possibleSol[k][nodeIdx] = true
                placedValue, delta2 = clusterValue(k, possibleSol, n, m, w_v, W_v, W)
                if placedValue <= B
                    tryCluster = false
                    goOn = false
                else
                    possibleSol[k][nodeToMove] = true
                    possibleSol[k][nodeIdx] = false
                end
            else
                possibleSol[clusterToTest][nodeToMove] = false
            end
            # try next cluster
            clusterToTest = mod(clusterToTest+1, K)
            if clusterToTest==0
                clusterToTest = K
            end
        end
        k += 1
    end
    return goOn
end

function affectRemainingNodes(fullSol::Bool, remainingNodes::Array{Int64}, sol, K::Int64, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64)::Bool
    success = true
    if fullSol==false
        for nodeIdx in remainingNodes
            println("TRYING TO PLACE NODE : ", nodeIdx)
            notPlaced = affectNode(nodeIdx, sol, K, B, n, m, w_v, W_v, W)
            if notPlaced==false
                println(nodeIdx, " placed")
            else
                success = false
            end
        end
    end
    return success
end

### TO GENERATE A SOLUTION BASED ON THE STATIC NODE WEIGHTS (DOES NOT ALWAYS SUCCEED)
function genSolByWeight(n::Int64, m::Int64, K::Int64, B::Int64, L::Int64, l, lh, w_v, W_v, W::Int64)
    println("dispatching nodes")
    fullSol, sol, remainingNodes = affectFirstNodes(n, m, K, B, L, l, lh, w_v, W_v, W)
    println("placing remaining nodes: ", remainingNodes)
    solFound = affectRemainingNodes(fullSol, remainingNodes, sol, K, B, n, m, w_v, W_v, W)
    println("SOLUTION: ", sol2Dto1D(sol))
    return solFound, sol
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

function heuristic(inputFile::String, timeLimit::Int64)
    include(inputFile)          # contains n, coordinates, lh[v], L, w_v, W_v, W, K, B,
    println("nb max de cluster : K = ", K)
    println("poids max d'un cluster : B = ", B)
    m = n*n
    l = generate_l(coordinates, n)

    start = time()

    println("LOOKING FOR FIRST SOLUTION")
    currentSol1D = genSolWithSolver(n, m, K, B, L, l, lh, w_v, W_v, W)
    currentSol2D = sol1Dto2D(currentSol1D, K, n)
    currentValue = partitionValue(currentSol1D, n, m, l, lh, L)
    println("found with value: ", currentValue)

    # getting the value, delta_2 and w_v2 for each cluster
    println("GETTING USEFUL DATA ABOUT INITIAL SOLUTION")
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

    # trying to move around vertices
    nodeIndices = [i for i in 1:n]
    shuffle!(nodeIndices)
    for nodeIdx in nodeIndices
        moved = simpleMove(nodeIdx, currentSol1D, currentSol2D, K, B, n, m, w_v, W_v, W, currentValue)
        if moved
            println("node ", nodeIdx, " changed of cluster")
        end
    end
    currentValue = partitionValue(currentSol1D, n, m, l, lh, L)
    println("value after search : ", currentValue)

    return currentSol1D, currentValue
end

    # generating a list of tuples of nodes with same w_v^2
    #couples = couplesWithSameW2(w_v2)
    #println("COUPLES : ", couples)

    #computationTime = time() - start

    #println("STARTING LOCAL SEARCH")
    #while computationTime < timeLimit
    #    # we try to switch two nodes from different clusters
    #    otherSol1D = switchTwoNodes(currentSol1D, couples)
    #    otherValue = partitionValue(otherSol1D, n, m, l, lh, L)
    #    if otherValue < currentValue
    #        currentValue = otherValue
    #        currentSol1D = otherSol1D
    #        println("better solution after ",  computationTime, "s with value : ", currentValue)
    #    end

    #    computationTime = time() - start
    #end

#    println("generated sol : ", currentSol1D)
#    return currentSol1D
#end

sol, val = heuristic("data/40_eil_6.tsp", 60)
println(sol)
