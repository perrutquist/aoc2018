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

function day3scan(s)
    m = match(r"#(\d*) @ (\d*),(\d*): (\d*)x(\d*)", s)
    map(x->parse(Int, x), m.captures)
end

function day3overlaps(v)
    mx = reduce((a,b) -> max.(a,b), v)
    M = zeros(Int, mx[2]+mx[4], mx[3]+mx[5])
    for c in v
        M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .+= 1
    end
    M
end

function day3part1(file)
    v = map(day3scan, readlines(file))
    sum(day3overlaps(v) .> 1)
end

function day3part2(file)
    v = map(day3scan, readlines(file))
    M = day3overlaps(v)
    for c in v
        all(M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .== 1) && return c[1]
    end
end

end # module

# Loop over puzzles and compute and print results
# (Reverse order, since we're most interested in the most recent result)
for d=25:-1:1, p=1:2
    f = Symbol("day$(d)part$(p)")
    if f ∈ names(AoC, all=true)
        r = @eval AoC.$f(string("input", $d, ".txt"))
        println("Day $d, part $p: ", r)
    end
end
