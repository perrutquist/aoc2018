module AoC

using InteractiveUtils # for copying to clipboard

include("solutions_per.jl") # Change this line to use your own code instead

"Scan one line of input for a given day"
scan(day::Integer, s) = scan(Val(Int(day)), s)
scan(::Val, s) = s # fallback

"Real all the input lines for a given day"
function getdata(day::Val{D}) where {D}
    file = "input$D.txt"
    isfile(file) || return nothing
    map(s->scan(day, s), eachline(file))
end
getdata(day::Integer) = getdata(Val(Int(day)))

"Solve the problem for a given day and part, given a vector of scanned lines."
solve(day::Integer, part::Integer, v) = solve(Val(Int(day)), Val(Int(part)), v)
solve(::Val, p, v) = nothing # fallback

const data = Array{Any,1}(undef, 25) # scanned input data
const solutions = Array{Any,2}(undef, 25, 2) # solutions

"""
   solveall()

Loop over puzzles, load data, compute and print results
"""
function solveall(;reverse=false, clipboard=false)
    t0 = time_ns()
    N = 0
    D = 0
    order = reverse ? Iterators.reverse : identity
    for d in order(1:25)
        global data[d] = getdata(d)
        if data[d] !== nothing
            D += 1
            println("\nDay $d: ", eltype(data[d]))
            for p in order(1:2)
                r = solve(d, p, data[d])
                global solutions[d, p] = r
                if r !== nothing
                    N += 1
                    print("Day $d, part $p: ")
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
solveall(reverse=true, clipboard=true)

end # module
