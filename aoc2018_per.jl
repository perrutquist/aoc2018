# Store the input data for day 1 as a file named `input1.txt`, and so on.
# Then run this file from the Julia REPL using the `include` command
#
# Re-running the `include` command will replace the entire module, re-scan
# all data, and re-run the computations. (No need to restart Julia.)
#
# At the REPL, each line of scanned data can be accessed as `AoC.data[day][line]`

module AoC

using InteractiveUtils

"Scan one line of input for a given day"
scan(day::Integer, s) = scanline(val(day), s)

"Real all the input lines for a given day"
function getdata(day::Val{D}) where {D}
    file = "input$D.txt"
    isfile(file) || return nothing
    map(s->scan(day, s), eachline(file))
end
getdata(day::Integer) = getdata(Val(day))

"Solve the problem for a given day and part, given a vector of scanned lines."
solve(day, part, v) = solve(Val(day), Val(part), v)

solve(::Val, ::Val, v) = nothing # fallback

# For each day we need a `scan` and two `solves` methods.
# day 0 can be copy-pasted as a starting point

# --------------------------------------------------------------------
day = 0

function scan(::Val{day}, s)
end

function solve(::Val{day}, ::Val{1}, v)
end

function solve(::Val{day}, ::Val{2}, v)
end

# --------------------------------------------------------------------
day = 1

scan(::Val{day}, s) = parse(Int, s)

function solve(::Val{day}, ::Val{1}, v)
    sum(v)
end

function solve(::Val{day}, ::Val{2}, v)
    f = 0
    s = Set([0])
    for d in Iterators.cycle(v)
        f += d
        f ∈ s && return f
        push!(s, f)
    end
end

# --------------------------------------------------------------------
day = 2

scan(::Val{day}, s) = Vector{Char}(s)

function solve(::Val{day}, ::Val{1}, v)
    match(v, n) = any(map(c -> sum(v .== c) == n, v))
    sum(match.(v, 2)) * sum(match.(v, 3))
end

function solve(::Val{day}, ::Val{2}, v)
    for a in v, b in v
        sum(a .!= b) == 1 && return String(a[a .== b])
    end
end

# --------------------------------------------------------------------
day = 3

function scan(::Val{day}, s)
    m = match(r"#(\d*) @ (\d*),(\d*): (\d*)x(\d*)", s)
    parse.((Int, Int, Int, Int, Int), Tuple(m.captures))
end

function day3overlaps(v)
    mx = reduce((a,b) -> max.(a,b), v)
    M = zeros(Int, mx[2]+mx[4], mx[3]+mx[5])
    for c in v
        M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .+= 1
    end
    M
end

function solve(::Val{day}, ::Val{1}, v)
    sum(day3overlaps(v) .> 1)
end

function solve(::Val{day}, ::Val{2}, v)
    M = day3overlaps(v)
    for c in v
        all(M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .== 1) && return c[1]
    end
end

# --------------------------------------------------------------------
# Loop over puzzles, load data, compute and print results
# (Reverse order, since we're usually interested in the most recent result)

const data = Array{Any,1}(undef, 25)

let t0 = time_ns(), N = 0, D = 0
    for d in 25:-1:1
        data[d] = getdata(d)
        if data[d] !== nothing
            D += 1
            println("\nDay $d: ", eltype(data[d]))
            for p in 2:-1:1
                r = solve(d, p, data[d])
                if r !== nothing
                    N += 1
                    print("Day $d, part $p: ")
                    if N == 1
                        InteractiveUtils.clipboard(string(r))
                        printstyled(r, bold=true)
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

end # module
