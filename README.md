# Advent of code in Julia (2018)

Store the input data for day 1 as a file named `input1.txt`, and so on.
Then run `include("aoc2018_per.jl")` from the Julia REPL.

Re-running the `include` command will replace the entire module, re-scan all data, and re-run the computations. (No need to restart Julia.)

At the REPL, each line of scanned data can be accessed as `AoC.data[day][line]`
