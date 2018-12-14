using OffsetArrays
using Plots

len = length # because i type "lenght" 50 % of the time when in a hurry.

# ----------------------------- Day  1 ---------------------------------------

scan(::Val{1}, line) = parse(Int, line)

function solve(day::Val{1}, ::Val{part}, lines) where part
    data = scan.(day, lines)

    if part ==1
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

function solve(day::Val{2}, ::Val{part}, lines) where part
    data = scan.(day, lines)

    if part == 1
        ismatch(v, n) = any(map(c -> sum(v .== c) == n, v))
        sum(ismatch.(data, 2)) * sum(ismatch.(data, 3))

    else # part 2
        for a in data, b in data
            sum(a .!= b) == 1 && return String(a[a .== b])
        end
    end
end

# ----------------------------- Day  3 ---------------------------------------

function scan(::Val{3}, line)::NTuple{5,Int}
    m = match(r"#(\d*) @ (\d*),(\d*): (\d*)x(\d*)", line)
    parse.((Int, Int, Int, Int, Int), Tuple(m.captures))
end

function solve(day::Val{3}, ::Val{part}, lines) where part
    data = scan.(day, lines)
    mx = reduce((a,b) -> max.(a,b), data)
    M = zeros(Int, mx[2]+mx[4], mx[3]+mx[5])
    for c in data
        M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .+= 1
    end

    if part == 1
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

function solve(day::Val{4}, ::Val{part}, lines) where part
    data = scan.(day, sort(lines)) # global D4 will be useless because it doesn't get sorted
    M = Dict{Int, Vector{Int}}()
    (t, id) = data[1](M, -1, 0)
    for f in data[2:end]
        (t, id) = f(M, t, id)
    end
    D = Dict([key => (part == 1 ? sum(val) : maximum(val)) for (key, val) in M])
    ix = argmax(D)
    ix*(argmax(M[ix])-1)
end

# ----------------------------- Day  5 ---------------------------------------

scan(::Val{5}, line) = Vector{Char}(line)

function solve(day::Val{5}, ::Val{part}, lines) where part
    polymer = scan(day, lines[1])

    function react(p, skip=nothing)
        o = Char[]
        for c in p
            uppercase(c) == skip && continue
            if !isempty(o) && abs(o[end] - c) == 32
                pop!(o)
            else
                push!(o, c)
            end
        end
        len(o)
    end

    if part == 1
        react(polymer)
    else # part 2
        minimum(react.((polymer,), 'A':'Z'))
    end
end

# ----------------------------- Day  6 ---------------------------------------

function scan(::Val{6}, line)
    m = match(r"(\d*), (\d*)", line)
    parse.((Int, Int), Tuple(m.captures))
end

function solve(day::Val{6}, ::Val{part}, lines; R = 10000) where part
    data = scan.(day, lines)
    manhattan(a, b) = sum(@. abs(a-b))
    (n,m) = reduce((a,b)->min.(a,b), data)
    (N,M) = reduce((a,b)->max.(a,b), data)
    D = OffsetArray{Int,2}(undef, n:N, m:M)
    C = Tuple.(CartesianIndices(D))

    if part == 1
        K = similar(D)
        dst = similar(D)
        D .= typemax(Int)
        for k in eachindex(data)
            p = data[k]
            @. dst = manhattan(C, (p,))
            @. K = ifelse(dst < D, k, ifelse(dst == D, -1, K))
            @. D = min(D, dst)
        end

        area(i) = any(K[[n, N], :] .== i) || any(K[:, [m, M]] .== i) ? -1 : sum(K .== i)
        maximum(area.(eachindex(data)))

    else # part 2
        # Note: The below code will not work for arbitrarily large R.
        D .= 0
        for p in data
            @. D += manhattan(C, (p,))
        end
        sum(D .< R)
    end
end

# ----------------------------- Day  7 ---------------------------------------

function scan(::Val{7}, line)::NTuple{2, Char}
    m = match(r"Step (.) must be finished before step (.) can begin.", line)
    first.(Tuple(m.captures))
end

function solve(day::Val{7}, ::Val{part}, lines; N = part==1 ? 1 : 5, D = 60) where part
    data = scan.(day, lines)

    o = Char[]
    l = unique(Iterators.flatten(data))
    second(x) = x[2]
    on = Vector{Union{Char, Nothing}}(nothing, N)
    til = fill(0, N)
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
                push!(o, on[i])
                til[i] = t + D + 1 + on[i] - 'A'
                filter!(!isequal(on[i]), l)
                isempty(l) && return part==1 ? String(o) : maximum(til)
            end
        end
    end
end

# ----------------------------- Day  8 ---------------------------------------

function scan(::Val{8}, line)
    parse.(Int, split(line, " "))
end

function solve(day::Val{8}, ::Val{part}, lines) where part
    data = scan(day, lines[1])

    function sm(i)
        n = data[i]
        m = data[i+1]
        s = 0
        r = i+2
        v = zeros(Int, n)
        for k in 1:n
            (r, sk, v[k]) = sm(r)
            s += sk
        end
        meta = data[r:r+m-1]
        (
        r+m,
        s+sum(meta),
        sum(k -> isempty(v) ? k : checkbounds(Bool, v, k) ? v[k] : 0, meta)
        )
    end
    sm(1)[part+1]
end

# ----------------------------- Day  9 ---------------------------------------

function scan(::Val{9}, line)::NTuple{2, Int}
    m = match(r"(\d*) players; last marble is worth (\d*) points", line)
    parse.(Int, Tuple(m.captures))
end

mutable struct Marble
    val::Int
    prev::Marble
    next::Marble
    function Marble(i)
        x = new(i)
        x.prev = x
        x.next = x
        x
    end
end

function solve(day::Val{9}, ::Val{part}, lines) where part
    (P, N) = scan(day, lines[1])

    function link(a,b)
        a.next = b
        b.prev = a
    end

    score = zeros(Int, P)
    c = Marble(0)
    for n in Marble.(1:(part==2 ? 100*N : N))
        if mod(n.val, 23)==0
            for k in 1:7
               c = c.prev
            end
            score[mod(n.val,P)+1] += c.val + n.val
            link(c.prev, c.next)
            c = c.next
        else
            c = c.next
            link(n, c.next)
            link(c, n)
            c = n
        end
    end
    maximum(score)
end

# ----------------------------- Day 10 ---------------------------------------

function scan(::Val{10}, line)::NTuple{4, Float64}
    m = match(r"position=< *(-?\d*), *(-?\d*)> velocity=< *(-?\d*), *(-?\d*)>", line)
    parse.(Float64, Tuple(m.captures))
end

function solve(day::Val{10}, ::Val{part}, lines) where part
    data = scan.(day, lines)

    p = getindex.(data, [1 2])
    v = getindex.(data, [3 4])

    x = [-v[:,2] ones(len(data))] \ p[:,2]

    t = round(x[1])
    @. p += t*v

    if part == 1
        plot(p[:,1], -p[:,2], seriestype=:scatter, aspect_ratio=1, legend=false)
        gui()
        return :solution_has_been_plotted
    else # part 2
        Int(t)
    end
end

# ----------------------------- Day 11 ---------------------------------------

function scan(::Val{11}, line)::Int
    parse.(Int, line)
end

function solve(day::Val{11}, ::Val{part}, lines) where part
    data = scan(day, lines[1])

    function fuel(x, y, d)
        r = x+10
        z = (r * y + d) * r
        mod(z ÷ 100, 10) - 5
    end

    D = fuel.(1:300, (1:300)', data)
    S1 = cumsum(D, dims=1)

    function bestsquare(r)
        S2 = S1[r:end, :]
        S2[2:end, :] .-= S1[1:end-r, :]
        S2 = cumsum(S2, dims=2)
        Q = S2[:, r:end]
        Q[:, 2:end] .-= S2[:, 1:end-r]
        findmax(Q)
    end

    if part == 1
        t=bestsquare(3)
        string(t[2][1], ",", t[2][2])

    else # part 2
        A = bestsquare.(1:300)
        i = argmax(first.(A))
        string(A[i][2][1], ",", A[i][2][2], ",", i)

    end
end

# ----------------------------- Day 12 ---------------------------------------

function solve(day::Val{12}, ::Val{part}, lines) where part
    G = part == 1 ? 20 : 1000

    s0 = Vector{Char}(lines[1][16:end]) .== '#'

    rules = Tuple(map(s -> s[10] == '#', sort(lines[2 .+ (1:32)], rev=true)))

    r(v) = rules[sum((16, 8, 4, 2, 1) .* v)+1]

    state = OffsetArray{Bool, 1}(undef, -2G:len(s0)+2G)
    state .= false
    state[0:len(s0)-1] .= s0
    s = copy(state)
    n = zeros(Int, G)
    for k in 1:G
        for j in firstindex(state)+2 : lastindex(state)-2
            state[j] = r((s[j-2], s[j-1], s[j], s[j+1], s[j+2]))
        end
        s .= state
        n[k] = sum( axes(state,1) .* state )
    end

    if part == 1
        return n[end]
    else
        # TODO: In the general case we'd need to find the period of the final cycle
        # For my input, it was 1.
        f(k) = n[G] + (n[G]-n[G-1])*(k-G)
        @assert f(G-10) == n[G-10]
        f(50000000000)
    end
end

# ----------------------------- Day 13 ---------------------------------------

function solve(day::Val{13}, ::Val{part}, lines) where part
    c = Vector{Char}.(lines)
    M = hcat(c...)
    M = permutedims(M,(2,1))
    M0 = copy(M)
    M0[M0 .== '<'] .= '-'
    M0[M0 .== '>'] .= '-'
    M0[M0 .== '^'] .= '|'
    M0[M0 .== 'v'] .= '|'
    T = zeros(Int,size(M))
    Mp = copy(M)

    while true
        Mp .= M
        for j in 1:size(M,2), i in 1:size(M,1)
            M[i,j] != Mp[i,j] && continue
            cp = M[i,j]
            t = T[i,j]
            if cp == '^'
                c = M0[i-1,j]
                if any(M[i-1,j] .== ('v', '<', '^', '>'))
                    part == 1 && return (j, i-1) .- (1, 1)
                    M[i,j] = M0[i,j]
                    M[i-1,j] = M0[i-1,j]
                    continue
                end
                (M[i-1,j], T[i-1,j]) = if c == '|'
                    ('^', t)
                elseif c == '/'
                    ('>', t)
                elseif c == '\\'
                    ('<', t)
                elseif c == '+'
                    if t == 0
                        ('<', 1)
                    elseif t == 1
                        ('^', 2)
                    elseif t == 2
                        ('>', 0)
                    end
                else
                    @show(i, j, c, t)
                    error("^")
                end
                M[i,j] = M0[i,j]
            elseif cp == '>'
                c = M0[i,j+1]
                if any(M[i,j+1] .== ('v', '<', '^', '>'))
                    part == 1 && return (j+1, i) .- (1, 1)
                    M[i,j] = M0[i,j]
                    M[i,j+1] = M0[i,j+1]
                    continue
                end
                (M[i,j+1], T[i,j+1]) = if c == '-'
                    ('>', t)
                elseif c == '/'
                    ('^', t)
                elseif c == '\\'
                    ('v', t)
                elseif c == '+'
                    if t == 0
                        ('^', 1)
                    elseif t == 1
                        ('>', 2)
                    elseif t == 2
                        ('v', 0)
                    end
                else
                    @show(i, j, c, t)
                    error(">")
                end
                M[i,j] = M0[i,j]
            elseif cp == 'v'
                c = M0[i+1,j]
                if any(M[i+1,j] .== ('v', '<', '^', '>'))
                    part == 1 && return (j, i+1) .- (1, 1)
                    M[i,j] = M0[i,j]
                    M[i+1,j] = M0[i+1,j]
                    continue
                end
                (M[i+1,j], T[i+1,j]) = if c == '|'
                    ('v', t)
                elseif c == '/'
                    ('<', t)
                elseif c == '\\'
                    ('>', t)
                elseif c == '+'
                    if t == 0
                        ('>', 1)
                    elseif t == 1
                        ('v', 2)
                    elseif t == 2
                        ('<', 0)
                    end
                else
                    @show(i, j, c, t)
                    error("^")
                end
                M[i,j] = M0[i,j]
            elseif cp == '<'
                c = M0[i,j-1]
                if any(M[i,j-1] .== ('v', '<', '^', '>'))
                    part == 1 && return (j-1, i) .- (1, 1)
                    M[i,j] = M0[i,j]
                    M[i,j-1] = M0[i,j-1]
                    continue
                end
                (M[i,j-1], T[i,j-1]) = if c == '-'
                    ('<', t)
                elseif c == '/'
                    ('v', t)
                elseif c == '\\'
                    ('^', t)
                elseif c == '+'
                    if t == 0
                        ('v', 1)
                    elseif t == 1
                        ('<', 2)
                    elseif t == 2
                        ('^', 0)
                    end
                else
                    @show(i, j, c, t)
                    error(">")
                end
                M[i,j] = M0[i,j]
            end #if
        end # for
        if part == 2
            s = 0
            for c in ('v', '<', '^', '>')
                s += sum(M .== c)
            end
            if s == 1
                for c in ('v', '<', '^', '>')
                    ix = findfirst( M .== c )
                    ix !== nothing && return string(ix[2]-1, ",", ix[1]-1)
                end
            end
        end

    end # while
end

# ----------------------------- Day 14 ---------------------------------------

function solve(day::Val{14}, ::Val{part}, lines) where part
    N = parse(Int, lines[1])
    v = Vector{Char}(lines[1]) .- '0'
    el = (1, 2)
    s = Int[3,7]
    while true
        r = getindex.((s,), el)
        n = +(r...)
        for c in (n >= 10 ? (1, n - 10) : (n,))
           push!(s, c)
           if part == 2 && len(s) >= len(v) && all(i -> s[end-len(v)+i] == v[i], 1:len(v))
               return len(s)-len(v)
           end
        end
        el = mod.(el .+ r, len(s)) .+ 1
        if part == 1 && len(s)>=N+10
            return string(s[N+1:N+10]...)
        end
    end
end
