module AoC

function day1part1(file)
    v = map(x->parse(Int, x), readlines(file))
    sum(v)
end

function day1part2(file)
    v = map(x->parse(Int, x), readlines(file))
    f = 0
    s = Set([0])
    for d in Iterators.cycle(v)
        f += d
        f ∈ s && return f
        push!(s, f)
    end
end

function day2part1(file)
    v = map(Vector{Char}, readlines(file))
    match(v, n) = any(map(c -> sum(v .== c) == n, v))
    sum(match.(v, 2)) * sum(match.(v, 3))
end

function day2part2(file)
    v = map(Vector{Char}, readlines(file))
    for a in v, b in v
        sum(a .!= b) == 1 && return String(a[a .== b])
    end
end

end # module

# Loop over puzzles and compute and print results
for d=1:24, p=1:2
    f = Symbol("day$(d)part$(p)")
    f ∈ names(AoC, all=true) || break
    r = @eval AoC.$f(string("input", $d, ".txt"))
    println("Day $d, part $p: ", r)
end
