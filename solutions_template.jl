# For each day, write a `scan` and a `solve` method.
# (Or skip the `scan` method and let `solve` start from a list of strings.)

# ----------------------------- Day  1 ---------------------------------------

scan(::Val{1}, line) = line

function solve(day::Val{1}, ::Val{part}, lines) where part
    data = scan.(day, lines)

    if part == 1

    else # part 2

    end
end

# ----------------------------- Day  2 ---------------------------------------

scan(::Val{2}, line) = line

function solve(day::Val{2}, ::Val{part}, lines) where part
    data = scan.(day, lines)

    if part == 1

    else # part 2

    end
end

# --------------------------------------------------------------------
# and so on...
