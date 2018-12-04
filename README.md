# Advent of code in Julia (2018)

This is a repository of my solutions to the AoC 2018 contest.

For maximum fun, try writing you own solutions before looking at somebody else's. However, there is some code in this repository that simply serves as a framework for looping over problems, loading data and printing results, and this code you might want to reuse.

## Usage

Copy `solutions_template.jl` and use as a starting point for writing your solutions.

Edit `aoc2018.jl` and modify the line that reads `include("solutions_per.jl")` (near the top) to instead use your solutions file.

For each day, two methods that need to be written:
* `scan(::Val{day}, s)` processes one line of input data and converts it to a more suitable format.
* `solve(::Val{day}, ::Val{part}, v)` takes a list `v` of whatever `scan` returned and computes the solution.

For days that are missing a `scan` method, `solve` will simply get a vector of the input strings.
If the `solve` method is missing or returns `nothing`, then no solution will be printed.

Store the input data for day `1` as a file named `input1.txt`, and so on. Then run `include("aoc2018_per.jl")` from the Julia REPL. The `solve` function will be solved for each day (1-25) and part (1-2).

Re-running the `include` command will replace the entire module, re-scan all data, and re-run the computations. (No need to restart Julia.)

At the REPL, each line of scanned data can be accessed as `AoC.data[day][line]`
