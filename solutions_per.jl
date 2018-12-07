using OffsetArrays
using StaticArrays

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
    ix = argmax(D)
    ix*(argmax(M[ix])-1)
end

# ----------------------------- Day  5 ---------------------------------------

scan(::Val{5}, line) = Vector{Char}(line)

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
    else # part 2
        same(a, c) = a == c || a == lowercase(c)
        f(c) = reactall!(polymer[@. !same(polymer,c)])
        minimum(f.('A':'Z'))
    end
end

# ----------------------------- Day  6 ---------------------------------------

function scan(::Val{6}, line)
    m = match(r"(\d*), (\d*)", line)
    parse.((Int, Int), Tuple(m.captures))
end

function solve(day::Val{6}, part, lines)
    data = scan.(day, lines)
    manhattan(a, b) = sum(@. abs(a-b))
    (n,m) = reduce((a,b)->min.(a,b), data)
    (N,M) = reduce((a,b)->max.(a,b), data)
    D = OffsetArray{Int,2}(undef, n:N, m:M)
    C = Tuple.(CartesianIndices(D))

    if part == Val(1)
        K = similar(D)
        dst = similar(D)
        D .= typemax(Int)
        for k in eachindex(data)
            p = data[k]
            @. dst = manhattan(C, (p,))
            @. K = ifelse(dst < D, k, ifelse(dst == D, -1, K))
            @. D = min(D, dst)
        end

        function area(i)
            any(K[n,:] .== i) && return -1
            any(K[N,:] .== i) && return -1
            any(K[:,m] .== i) && return -1
            any(K[:,M] .== i) && return -1
            sum(K .== i)
        end

        maximum(area.(eachindex(data)))

    else # part 2
        R = 10000 # sum of distances must be less than this
        # Note: The below code will not work for arbitrarily large R.

        D .= 0
        for p in data
            @. D += manhattan(C, (p,))
        end
        sum(D .< R)
    end
end

# ----------------------------- Day  7 ---------------------------------------

function scan(::Val{7}, line)
    m = match(r"Step (.) must be finished before step (.) can begin.", line)
    first.(Tuple(m.captures))
end

function solve(day::Val{7}, part, lines)
    data = scan.(day, lines)
    o = ""
    l = unique(Iterators.flatten(data))
    second(x) = x[2]

    if part === Val(1)
        while true
            l2 = unique(second.(data))
            a = minimum(setdiff(l, l2))
            o *= a
            filter!(d -> first(d) != a, data)
            filter!(!isequal(a), l)
            isempty(l) && return o
        end

    else # part 2
        on = Vector{Union{Char, Nothing}}(nothing, 5)
        til = fill(0, 5)
        while true
            t = minimum(til)
            for i in findall(til .== t)
                a = on[i]
                filter!(d -> first(d) != a, data)
                l2 = unique(second.(data))
                s = setdiff(l, l2)
                if isempty(s)
                    on[i] = nothing
                    til[i] = minimum(til[on.!=nothing])
                else
                    on[i] = minimum(s)
                    til[i] = t + 61 + on[i]-'A'
                    filter!(!isequal(on[i]), l)
                    isempty(l) && return maximum(til)
                end
            end
        end
    end
end
