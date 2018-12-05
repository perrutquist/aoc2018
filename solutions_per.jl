# ----------------------------- Day  1 ---------------------------------------

scan(::Val{1}, line) = parse(Int, line)

function solve(day::Val{1}, part, lines)
    data = scan.(day, lines)
    if part === Val(1)
        sum(data)
    else
        f = 0
        s = Set([0])
        for d in Iterators.cycle(data)
            f += d
            f ∈ s && return f
            push!(s, f)
        end
    end
end

# ----------------------------- Day  2 ---------------------------------------

scan(::Val{2}, line) = Vector{Char}(line)

function solve(day::Val{2}, part, lines)
    data = scan.(day, lines)
    if part === Val(1)
        ismatch(v, n) = any(map(c -> sum(v .== c) == n, v))
        sum(ismatch.(data, 2)) * sum(ismatch.(data, 3))
    else # part 2
        for a in data, b in data
            sum(a .!= b) == 1 && return String(a[a .== b])
        end
    end
end

# ----------------------------- Day  3 ---------------------------------------

function scan(::Val{3}, line)
    m = match(r"#(\d*) @ (\d*),(\d*): (\d*)x(\d*)", line)
    parse.((Int, Int, Int, Int, Int), Tuple(m.captures))
end

function solve(day::Val{3}, part, lines)
    data = scan.(day, lines)
    mx = reduce((a,b) -> max.(a,b), data)
    M = zeros(Int, mx[2]+mx[4], mx[3]+mx[5])
    for c in data
        M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .+= 1
    end
    if part === Val(1)
        sum(M .> 1)
    else # part 2
        for c in data
            all(M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .== 1) && return c[1]
        end
    end
end

# ----------------------------- Day  4 ---------------------------------------

function scan(::Val{4}, line)
    m = match(r"\[\d*-\d*-\d* \d*:(\d*)\] (.*)", line)
    s = parse(Int, m.captures[1])

    if (g = match(r"Guard #(\d*) begins shift", m.captures[2])) !== nothing
        id = parse(Int, g.captures[1])
        (M, t, xx)->begin
            @assert t==-1
            haskey(M,id) || (M[id] = zeros(Int, 60))
            (-1, id)
        end
    elseif m.captures[2]=="falls asleep"
        (M, t, id)->(s, id)
    else # wakes up
        (M, t, id)->(M[id][(t+1):s] .+= 1; (-1, id))
    end
end

function solve(day::Val{4}, part, lines)
    data = scan.(day, sort(lines)) # global D4 will be useless because it doesn't get sorted
    M = Dict{Int, Vector{Int}}()
    (t, id) = data[1](M, -1, 0)
    for f in data[2:end]
        (t, id) = f(M, t, id)
    end
    D = Dict([key => (part === Val(1) ? sum(val) : maximum(val)) for (key, val) in M])
    (xx, ix) = findmax(D)
    (xx, mn) = findmax(M[ix])
    ix*(mn-1)
end

# ----------------------------- Day  5 ---------------------------------------

scan(::Val{5}, line) = Vector{UInt8}(line)

function solve(day::Val{5}, part, lines)
    polymer = scan.(day, lines)[1]

    function reactall!(p)
        r(i, m) = i%2 == m && abs( p[i] - p[i+1] ) == 32
        reacting(i, m) = (i<length(p) && r(i, m)) || (i > 1 && r(i-1, m))
        while true
            l = length(p)
            for m in 0:1
                deleteat!(p, reacting.(1:length(p), m))
            end
            length(p) == l && return l
        end
    end

    if part === Val(1)
        reactall!(polymer)
    else
        same(a, c) = a == c || a == lowercase(c)
        f(c) = reactall!(polymer[@. !same(polymer,c)])
        minimum(f.('A':'Z'))
    end
end
