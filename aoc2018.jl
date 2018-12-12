module AoC

using InteractiveUtils # for copying to clipboard
using Printf
using BenchmarkTools

include("solutions_per.jl") # Change this line to use your own code instead

"""
    main(a::Symbol, v)

Assign `v` to a variable in the `Main` module, with name given by the symbol `a`,
and display it to stdout.

This can be useful for debugging, but should probably never be used in
production code.
"""
function main(a::Symbol, v, verbose)
    Main.eval( :( $a = $v ) )
    if verbose
        println(a, " =")
        display(v)
    end
end

"""
    @main x

Copy a variable (e.g. `x`) into the `Main` module, and display it to stdout.

This can be useful for debugging, but should probably never be used in
production code.
"""
macro main(a::Symbol, verbose::Bool=true)
    quote
        main($(QuoteNode(a)), $(esc(a)), $verbose)
    end
end

"Parse one line of input for a given day and convert to a more useful form"
scan(day::Integer, s) = scan(Val(Int(day)), s)
scan(::Val, s) = s # fallback

"Read all the input lines for a given day into a list"
function getlines(day)
    file = "input$day.txt"
    isfile(file) || return nothing
    readlines(file)
end

"Solve the problem for a given day and part, given a vector of scanned lines."
solve(day::Integer, part::Integer, v; kwargs...) = solve(Val(Int(day)), Val(Int(part)), v; kwargs...)
solve(::Val, p, v; kwargs...) = nothing # fallback

const lines = Array{Any,1}(undef, 25) # raw input as list of strings
const solution = Array{Any,2}(undef, 25, 2)
const timing = zeros(25, 2)

"""
   solveall()

Loop over puzzles, load data, compute and print results
"""
function solveall(;reverse=false, clipboard=false, setglobals=false, twice=false)
    t0 = time_ns()
    N = 0
    D = 0
    order = reverse ? Iterators.reverse : identity
    for day in order(1:25)
        ld = getlines(day)
        if ld !== nothing
            lines[day] = ld
            D += 1
            if setglobals
                Dn = Symbol("D$day")
                dd = scan.(day, deepcopy(ld))
                println("\nDay $day: ", eltype(dd))
                Main.eval(:( $Dn = $dd ))
            end
            for part in order(1:2)
                td0 = time_ns()
                r = solve(day, part, deepcopy(ld))
                global solution[day, part] = r
                t = time_ns() - td0 # nanoseconds
                global timing[day, part] = t*1e-9 # seconds
                if r !== nothing
                    N += 1
                    print("Day $day, part $part: ")
                    if N == 1 && clipboard
                        printstyled(r, bold=true)
                        InteractiveUtils.clipboard(string(r))
                        print("  (copied to clipboard)")
                    else
                       print(r)
                    end
                    Printf.@printf "  Time: %.3f ms" t*1e-6

                    if twice
                        td1 = time_ns()
                        r2 = solve(day, part, deepcopy(ld))
                        t2 = time_ns() - td1
                        @assert r == r2
                        Printf.@printf ", 2:nd time: %.3f ms." t2*1e-6
                    end

                    println()
                end
            end
        end
    end
    if D != 0
        Printf.@printf("\nSolved %d problems in %.6f seconds.\n", N, sum(timing))
        println("First timings may include JIT-compilation time, 2:nd timings measure only execution.")
        #Printf.@printf("Total time, including file reads: %.2f seconds.\n", 1e-9*(time_ns()-t0))
    else
        println("No data files found!")
    end
end

"""
   benchmarks()

Loop over puzzles and print benchmark results.
Each benchmark is calculated by taking the minimum time over a large number of runs.
"""
function benchmarks()
    for day in 1:25
        ld = getlines(day)
        if ld !== nothing
            lines[day] = ld
            for part in 1:2
                r = solve(day, part, deepcopy(ld))
                if r !== nothing
                    print("Day $day, part $part: ")
                    @btime solve($day, $part, l) setup=(l=deepcopy($ld))
                end
            end
        end
    end
end

# Solve all puzzles.
# Reverse order, since we're usually interested in the most recent result.
# Copy to clipboard to make it easier to submit to contest.
solveall(reverse=true, clipboard=true, setglobals=true, twice=true)

end # module
