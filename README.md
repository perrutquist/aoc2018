# Advent of Code (2018) in Julia

This is a repository of my solutions to the AoC 2018 contest.

For maximum fun, try writing you own solutions before looking at somebody else's. However, there is some code in this repository that simply serves as a framework for looping over problems, loading data and printing results, and this code you might want to reuse.

## Usage

Copy `solutions_template.jl` and use as a starting point for writing your solutions.

Edit `aoc2018.jl` and modify the line that reads `include("solutions_per.jl")` (near the top) to instead use your solutions file.

For each day, a `solve` method needs to be written.
```
solve(::Val{day}, ::Val{part}, list)
```
This function takes a list of strings (input lines) and computes the solution.

The `solve` methods for each day will usually consist of some code that is common to the two parts, followed by an `if` statement that computes the solution for each part. However, it is also possible to write completely separate methods for the two parts.

It is convenient to start each `solve` method with `data = scan.(day, lines)`, where `scan` is a function that converts each line into a more useful format. For example:
```
scan(::Val{1}, line) = parse(Int, line)
```

Store the input data for day `1` as a file named `input1.txt`, and so on. Then run `include("aoc2018.jl")` from the Julia REPL. The `solve` function will be solved for each day (1-25) and part (1-2).

Re-running the `include` command will replace the entire module, re-scan all data, and re-run the computations. (No need to restart Julia.)

The global variables `D1`, `D2`... will contain the input data (as returned by `scan`).
