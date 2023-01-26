include("dualisation.jl")
include("coupes.jl")

using CSV, DataFrames

function benchmarkDualisation()
    # cd("./data/")
    files = readdir("./data_dual")
    times = []
    values = []
    solutions = []
    for file in files
        println("Processing " * file)
        xStar, clusters, value, compTime = master_pb("data_dual/" * file)
        push!(times, compTime)
        push!(values, value)
        push!(solutions, clusters)
    end
    df = DataFrame(File = files,
                   Time = times,
                   Values = values,
                   Solutions = solutions)
    CSV.write("benchmarkDualisation.csv", df, delim=";")
end

function benchmarkCoupes(timeLimit::Int64)
    # cd("./data/")
    files = readdir("./data_cuts")
    statuses = []
    times = []
    values = []
    solutions = []
    for file in files
        println("Processing " * file)
        clusters, value, lower_bound, compTime, status = solveByCuts("data_cuts/" * file, timeLimit)
        push!(statuses, status)
        push!(times, compTime)
        push!(values, value)
        push!(solutions, clusters)
    end
    df = DataFrame(File = files,
                   Time = times,
                   Values = values,
                   Solutions = solutions)
    CSV.write("benchmarkCoupes.csv", df, delim=";")
end

# benchmarkDualisation()
