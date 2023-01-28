function removeCouplesInSameCluster(vectorSol, couples)
    couplesToSwich = []
    for c in couples
        if vectorSol[c[1]] != vectorSol[c[2]]
            push!(couplesToSwich, c)
        end
    end
    return couplesToSwich
end

function switchTwoNodes(vectorSol, couples)
    couplesToSwich = removeCouplesInSameCluster(vectorSol, couples)
    randomIdx = rand(1:length(couplesToSwich))
    c = couplesToSwich[randomIdx]
    vectorSol[c[1]] = vectorSol[c[2]]
    vectorSol[c[2]] = vectorSol[c[1]]
    return vectorSol
end