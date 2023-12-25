include("src/MetaSchwinger.jl")
using .MetaSchwinger

if length(ARGS) == 0
    error("""
    Use input file:
    Like,
    julia run.jl parameters.jl
    """)
end

function runtest()
    run_sim(ARGS[1])
end

@time runtest()