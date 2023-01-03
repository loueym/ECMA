function distance(Xp1::Float64, Yp1::Float64, Xp2::Float64, Yp2::Float64)::Float64
    return sqrt((Xp1-Xp2)*(Xp1-Xp2) + (Yp1-Yp2)*(Yp1-Yp2))
end

function generate_l(coordinates::Array{Float64,2}, n::Int)::Vector{Float64}
    l = Vector{Float64}()
    for v1 in 1:n
        for v2 in 1:n
            append!(l, distance(coordinates[v1,1], coordinates[v1,2], coordinates[v2,1], coordinates[v2,2]))
        end
    end
    return l
end

function node1(e::Int64, n::Int64)::Int64
    return (e-1)÷n+1
end

function node2(e::Int64, n::Int64)::Int64
    return (e-1)%n+1
end
