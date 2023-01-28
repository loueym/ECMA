include("dualisation.jl")
include("plansCoupants.jl")
include("branchAndCut.jl")

using CSV, DataFrames

function benchmarkDualisation(timeLimit::Int)
    # cd("./data/")
    files = readdir("./data_bnc")
    gaps = []
    times = []
    values = []
    solutions = []
    for file in files
        println("Processing " * file)
        xStar, clusters, value, compTime, gap = solveByDualisation("data_bnc/" * file, timeLimit)
        push!(gaps, gap)
        push!(times, compTime)
        push!(values, value)
        push!(solutions, clusters)
    end
    df = DataFrame(File = files,
                   Time = times,
                   Values = values,
                   Gaps = gaps,
                   Solutions = solutions)
    CSV.write("benchmarkDualisation.csv", df, delim=";", append=true)
end

function benchmarkBnC(timeLimit::Int64)
    # cd("./data/")
    files = readdir("./data_bnc")
    statuses = []
    gaps = []
    times = []
    values = []
    solutions = []
    for file in files
        println("Processing " * file)
        clusters, t_star, value, compTime, gap, status = solveByBnC("data_bnc/" * file, timeLimit)
        push!(statuses, status)
        push!(gaps, gap)
        push!(times, compTime)
        push!(values, value)
        push!(solutions, clusters)
    end
    df = DataFrame(File = files,
                   Time = times,
                   Values = values,
                   Gaps = gaps,
                   Solutions = solutions)
    CSV.write("benchmarkBnC.csv", df, delim=";", append=true)
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
