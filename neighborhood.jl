include("common.jl")

### TO COMPUTE THE VALUE OF THE GIVEN PARTITION
function partitionValue(sol1D, n::Int64, m::Int64, l, lh, L::Int64)::Float64
    x = xFromPartition(sol1D, n, m)

    solVal = Model(CPLEX.Optimizer)
    @variable(solVal, delta1[e in 1:m] >= 0)

    @constraint(solVal, sum(delta1[e] for e in 1:m) <= L)
    @constraint(solVal, [e in 1:m], delta1[e] <= 3)

    @objective(solVal, Max, sum(x[e]*(1/2*l[e]+delta1[e]*(lh[node1(e, n)]+lh[node2(e, n)])) for e in 1:m))

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(solVal, "CPX_PARAM_SCRIND", 0)
    optimize!(solVal)

    feasiblefound = primal_status(solVal) == MOI.FEASIBLE_POINT
    if feasiblefound
        value = JuMP.objective_value(solVal)
    end

    return value
end

### TO COMPUTE THE VALUE OF CLUSTER k IN THE GIVEN PARTITION
function clusterValue(k::Int64, partition2D, n::Int64, m::Int64, w_v, W_v, W::Int64)
    clusVal = Model(CPLEX.Optimizer)
    @variable(clusVal, delta2[v in 1:n] >= 0)

    @constraint(clusVal, sum(delta2[v] for v in 1:n) <= W)
    @constraint(clusVal, [v in 1:n], delta2[v] <= W_v[v])

    @objective(clusVal, Max, sum(w_v[v]*partition2D[k][v]*(1+delta2[v]) for v in 1:n))

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(clusVal, "CPX_PARAM_SCRIND", 0)
    optimize!(clusVal)

    feasiblefound = primal_status(clusVal) == MOI.FEASIBLE_POINT
    if feasiblefound
        delta2Star = JuMP.value.(delta2)
        clusterVal = JuMP.objective_value(clusVal)
    end

    return clusterVal, delta2Star
end

# returns a vector with couples that are not in the same cluster
function removeSameClusterCouples(vectorSol::Vector{Int64}, couples)
    couplesToSwitch = Vector()
    for c in couples
        if vectorSol[c[1]] != vectorSol[c[2]]
            push!(couplesToSwitch, c)
        end
    end
    return couplesToSwitch
end

# swaps a couple at random, making sure both nodes are not in the same cluster
function switchTwoNodes(sol1D::Vector{Int64}, sol2D, currentValue::Float64, couples, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64, l, lh, L)
    # couples is a vector of tuples of nodes that are not in the same cluster and that could be swiched
    shuffle!(couples)
    nodesChanged = Vector{Int64}()
    for i in 1:length(couples)
        c = couples[i]
        n1, n2 = c[1], c[2]
        if !(n1 in nodesChanged) && !(n2 in nodesChanged)
            k1, k2 = sol1D[n1], sol1D[n2]
            sol2D[k2][n2] = false
            sol2D[k1][n2] = true
            sol2D[k1][n1] = false
            sol2D[k2][n1] = true
            cluster1Val, delta21 = clusterValue(k1, sol2D, n, m, w_v, W_v, W)
            cluster2Val, delta22 = clusterValue(k2, sol2D, n, m, w_v, W_v, W)
            if cluster1Val <= B && cluster2Val <= B
                sol1D[n1], sol1D[n2] = k2, k1
                movedValue = partitionValue(sol1D, n, m, l, lh, L)
                if movedValue < currentValue
                    # we found a better solution
                    println("changing sol ! current = ", currentValue, " and new = ", movedValue)
                    append!(nodesChanged, n1)
                    append!(nodesChanged, n2)
                    currentValue = movedValue
                else
                    # undo changes
                    sol1D[n1], sol1D[n2] = k1, k2
                    sol2D[k2][n2] = true
                    sol2D[k1][n2] = false
                    sol2D[k1][n1] = true
                    sol2D[k2][n1] = false
                end
            else
                # undo changes
                sol2D[k2][n2] = true
                sol2D[k1][n2] = false
                sol2D[k1][n1] = true
                sol2D[k2][n1] = false
            end
        end
    end
    println("current value is ", currentValue, " and partition value is ", partitionValue(sol1D, n, m, l, lh, L))
end

# moves node nodeIdx to another cluster if it enhances the solution
function simpleMove(nodeIdx::Int64, sol1D::Array{Int64}, sol2D, K::Int64, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64, l, lh, L::Int64, currentValue::Float64)::Bool
    originCluster = sol1D[nodeIdx]
    destinationCluster = 0
    bestValue = currentValue
    # k is the cluster to put the node in
    if originCluster==1
        k = 2
    else
        k = 1
    end
    while k<=K
        # test simple move
        sol2D[k][nodeIdx] = true
        sol1D[nodeIdx] = k
        movedClusterValue, delta2 = clusterValue(k, sol2D, n, m, w_v, W_v, W)
        if movedClusterValue <= B
            movedValue = partitionValue(sol1D, n, m, l, lh, L)
            if movedValue < bestValue
                destinationCluster = k
                bestValue = movedValue
            end
        end
        # get back to original partition
        sol2D[k][nodeIdx] = false
        # test the next cluster
        k += 1
        if k==originCluster
            k+=1
        end
    end
    if destinationCluster!=0
        sol2D[originCluster][nodeIdx] = false
        sol2D[destinationCluster][nodeIdx] = true
        sol1D[nodeIdx] = destinationCluster
        return true
    else
        sol1D[nodeIdx] = originCluster
        return false
    end
end

function nodesWithSameW(w_v, allowedGap::Float64)
    nodesBywv = Dict{}()
    for i in 1:length(w_v)
        w_v_i = w_v[i]
        if w_v_i in keys(nodesBywv)
            append!(nodesBywv[w_v_i], i)
        else
            nodesBywv[w_v_i] = Vector{Int64}([i])
        end
        for w_v_val in keys(nodesBywv)
            if (w_v_val != w_v_i) && (abs(w_v_i - w_v_val)/w_v_val < allowedGap)
                append!(nodesBywv[w_v_val], i)
            end
        end
    end
    return nodesBywv
end

function changeSol2D(nodes::Vector{Int64}, allK::Vector{Int64}, sol2D, undo::Bool)
    for i in 1:length(nodes)
        ni = nodes[i]
        k = allK[i]
        if i < length(nodes)
            nextK = allK[i+1]
        else
            nextK = allK[1]
        end
        sol2D[k][ni] = undo
        sol2D[nextK][ni] = !(undo)
    end
end

function changeSol1D(nodes::Vector{Int64}, allK::Vector{Int64}, sol1D, undo::Bool)
    for i in 1:length(nodes)
        ni = nodes[i]
        if undo
            sol1D[ni] = allK[i]
        else
            if i < length(nodes)
                nextK = allK[i+1]
            else
                nextK = allK[1]
            end
            sol1D[ni] = nextK
        end
    end
end

function swichNodes(sol1D::Vector{Int64}, sol2D, currentValue::Float64, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64, l, lh, L::Int64)
    for gapValue in [0, 0.05, 0.1, 0.15]
        nodesBywv = nodesWithSameW(w_v, gapValue)
        nodesBywvValues = shuffle!(collect(values(nodesBywv)))
        for nodes in nodesBywvValues
            shuffle!(nodes)
            nbNodes = length(nodes)
            if nbNodes > 1
                allK = Vector{Int64}([sol1D[i] for i in nodes])
                changeSol2D(nodes, allK, sol2D, false)
                
                allLessThanB = true
                for i in 1:nbNodes
                    clusterVal, delta2 = clusterValue(allK[i], sol2D, n, m, w_v, W_v, W)
                    if clusterVal > B
                        allLessThanB = false
                        break
                    end
                end
                if allLessThanB
                    changeSol1D(nodes, allK, sol1D, false)
                    movedValue = partitionValue(sol1D, n, m, l, lh, L)
                    if movedValue < currentValue
                        # we found a better solution
                        println("changing sol ! current = ", currentValue, " and new = ", movedValue)
                        currentValue = movedValue
                    else
                        # undo changes
                        changeSol1D(nodes, allK, sol1D, true)
                        changeSol2D(nodes, allK, sol2D, true)
                    end
                else
                    # undo changes
                    changeSol2D(nodes, allK, sol2D, true)
                end
            end
        end
    end
end