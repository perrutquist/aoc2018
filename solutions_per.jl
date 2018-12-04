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

function solve(::Val{day}, part, v)
    mx = reduce((a,b) -> max.(a,b), v)
    M = zeros(Int, mx[2]+mx[4], mx[3]+mx[5])
    for c in v
        M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .+= 1
    end
    if part === Val(1)
        return sum(M .> 1)
    else # part 2
        for c in v
            all(M[c[2]+1:c[2]+c[4], c[3]+1:c[3]+c[5]] .== 1) && return c[1]
        end
    end
end

# --------------------------------------------------------------------
day = 4

function scan(::Val{day}, s)
    m = match(r"\[(\d*)-(\d*)-(\d*) (\d*):(\d*)\] (.*)", s)
    s = parse.(Int, m.captures[1:5])
    g = match(r"Guard #(\d*) begins shift", m.captures[6])
    id = g === nothing ? 0 : parse(Int, g.captures[1])
    (366*s[1] + 31*s[2] + s[3]+(s[4]==23), s[5]-60*(s[4]==23), id, m.captures[6]=="falls asleep")
end

function solve(::Val{day}, part, v)
    sort!(v, by=i->256*i[1]+i[2])
    id = 0
    t = 0
    N = maximum(x->x[3], v)
    s = zeros(Int, N, 60)
    sleeps = false
    for (xx, m, i, sl) in v
        if i > 0
            if sleeps # last guard slept at end of hour
                # This never happened in the input data
                s[id, t+1:end] .+= 1
            end
            id = i
            sleeps = false
        elseif sl
            t = m
            sl = true
        else # wakes up
            s[id, t+1:m] .+= 1
            sleeps = false
        end
    end
    if part == Val(1)
        (xx, ix) = findmax(vec(sum(s, dims=2)))
        (xx, mn) = findmax(s[ix,:])
        return (ix*(mn-1))
    else # part 2
        (xx, ix) = findmax(s)
        return ix[1]*(ix[2]-1)
    end
end
