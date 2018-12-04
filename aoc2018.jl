module AoC

using InteractiveUtils # for copying to clipboard

include("solutions_per.jl") # Change this line to use your own code instead

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
solve(day::Integer, part::Integer, v) = solve(Val(Int(day)), Val(Int(part)), v)
solve(::Val, p, v) = nothing # fallback

const lines = Array{Any,1}(undef, 25) # raw input as list of strings
const solutions = Array{Any,2}(undef, 25, 2) # solutions

"""
   solveall()

Loop over puzzles, load data, compute and print results
"""
function solveall(;reverse=false, clipboard=false, setglobals=false)
    t0 = time_ns()
    N = 0
    D = 0
    order = reverse ? Iterators.reverse : identity
    for day in order(1:25)
        ld = getlines(day)
        if ld !== nothing
            D += 1
            if setglobals
                Dn = Symbol("D$day")
                dd = scan.(day, ld)
                println("\nDay $day: ", eltype(dd))
                Main.eval(:( $Dn = $dd ))
            end
            for part in order(1:2)
                r = solve(day, part, ld)
                global solutions[day, part] = r
                if r !== nothing
                    N += 1
                    print("Day $day, part $part: ")
                    if N == 1 && clipboard
                        printstyled(r, bold=true)
                        InteractiveUtils.clipboard(string(r))
                        println(" (copied to clipboard)")
                    else
                       println(r)
                    end
                end
            end
        end
    end
    if D != 0
        println("\nSolved $N problems in $(1e-9*(time_ns()-t0)) seconds (including JIT-compile-time).")
    else
        println("No data files found!")
    end
end

# Solve all puzzles.
# Reverse order, since we're usually interested in the most recent result.
# Copy to clipboard to make it easier to submit to contest.
solveall(reverse=true, clipboard=true, setglobals=true)

end # module
