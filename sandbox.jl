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
