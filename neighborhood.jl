include("common.jl")

### TO COMPUTE THE VALUE OF THE GIVEN PARTITION
function partitionValue(vectorPartition, n::Int64, m::Int64, l, lh, L::Int64)::Float64
    x = xFromPartition(vectorPartition, n, m)

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
function switchTwoNodes(vectorSol::Vector{Int64}, couples)
    couplesToSwitch = removeSameClusterCouples(vectorSol, couples)
    randomIdx = rand(1:length(couplesToSwitch))
    c = couplesToSwitch[randomIdx]
    swappedVectorSol = copy(vectorSol)
    swappedVectorSol[c[1]] = vectorSol[c[2]]
    swappedVectorSol[c[2]] = vectorSol[c[1]]
    return swappedVectorSol
end

# moves node nodeIdx to another cluster if it enhances the solution
function simpleMove(nodeIdx::Int64, sol1D::Array{Int64}, sol2D, K::Int64, B::Int64, n::Int64, m::Int64, w_v, W_v, W::Int64, currentValue::Float64)
    originCluster = sol1D[nodeIdx]
    # k is the cluster to put the node in
    if originCluster==1
        k = 2
    else
        k = 1
    end
    goOn = true
    while k<=K && goOn
        # test simple move
        sol2D[k][nodeIdx] = true
        movedClusterValue, delta2 = clusterValue(k, sol2D, n, m, w_v, W_v, W)
        if movedClusterValue <= B
            sol1D[nodeIdx] = k
            movedValue = partitionValue(sol1D, n, m, l, lh, L)
            if movedValue < currentValue
                goOn = false
                sol2D[originCluster][nodeIdx] = false
            else
                sol1D[nodeIdx] = originCluster
                sol2D[k][nodeIdx] = false
                sol2D[originCluster][nodeIdx] = true
            end
        else
            sol2D[k][nodeIdx] = false
        end
        k += 1
        if k==originCluster
            k+=1
        end
    end
    return !goOn
end
