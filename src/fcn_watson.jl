using DrWatson

"""
    resultsdir(line)

Path to results directory.
"""
resultsdir(args...) = projectdir("results", args...)


"""
    benchdir(line)

Path to benchmark results directory.
"""
benchdir(args...) = projectdir("results", "bench", args...)


"""
    modelsdir(line)

Path to models directory.
"""
modelsdir(args...) = projectdir("src", "models", args...)


"""
    flushln(line)

Print line and flush standard output.

# Arguments
- `line::String`: string to flush to standard output
"""
function flushln(line::String)
    println(line)
    flush(stdout)
end
