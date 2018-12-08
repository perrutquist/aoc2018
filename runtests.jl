# TODO: Set up a proper package

# This script assumes that the AoC module has been loaded

using Test

lines(str) = chomp.(split(str, "\n")[1:end-1])

L = Vector{Vector{String}}(undef, 25)

day = 1

L[day] = lines("""
+1
-2
+3
+1
""")

@test AoC.solve(day, 1, L[day]) == 3
@test AoC.solve(day, 2, L[day]) == 2


day = 2

L[day] = lines("""
abcdef
bababc
abbcde
abcccd
aabcdd
abcdee
ababab
""")

@test AoC.solve(day, 1, L[day]) == 12

L[day] = lines("""
abcde
fghij
klmno
pqrst
fguij
axcye
wvxyz
""")

@test AoC.solve(day, 2, L[day]) == "fgij"


day = 3

L[day] = lines("""
#1 @ 1,3: 4x4
#2 @ 3,1: 4x4
#3 @ 5,5: 2x2
""")

@test AoC.solve(day, 1, L[day]) == 4
@test AoC.solve(day, 2, L[day]) == 3


day = 4

L[day] = lines("""
[1518-11-01 00:00] Guard #10 begins shift
[1518-11-01 00:05] falls asleep
[1518-11-01 00:25] wakes up
[1518-11-01 00:30] falls asleep
[1518-11-01 00:55] wakes up
[1518-11-01 23:58] Guard #99 begins shift
[1518-11-02 00:40] falls asleep
[1518-11-02 00:50] wakes up
[1518-11-03 00:05] Guard #10 begins shift
[1518-11-03 00:24] falls asleep
[1518-11-03 00:29] wakes up
[1518-11-04 00:02] Guard #99 begins shift
[1518-11-04 00:36] falls asleep
[1518-11-04 00:46] wakes up
[1518-11-05 00:03] Guard #99 begins shift
[1518-11-05 00:45] falls asleep
[1518-11-05 00:55] wakes up
""")

@test AoC.solve(day, 1, L[day]) == 240
@test AoC.solve(day, 2, L[day]) == 4455


day = 5

L[day] = lines("""
dabAcCaCBAcCcaDA
""")

@test AoC.solve(day, 1, L[day]) == 10
@test AoC.solve(day, 2, L[day]) == 4


day = 6

L[day] = lines("""
1, 1
1, 6
8, 3
3, 4
5, 5
8, 9
""")

@test AoC.solve(day, 1, L[day]) == 17
@test AoC.solve(day, 2, L[day]; R=32) == 16


day = 7

L[day] = lines("""
Step C must be finished before step A can begin.
Step C must be finished before step F can begin.
Step A must be finished before step B can begin.
Step A must be finished before step D can begin.
Step B must be finished before step E can begin.
Step D must be finished before step E can begin.
Step F must be finished before step E can begin.
""")

@test AoC.solve(day, 1, L[day]) == "CABDFE"
@test AoC.solve(day, 2, L[day], N=2, D=0) == 15


day = 8

L[day] = lines("""
2 3 0 3 10 11 12 1 1 0 1 99 2 1 1 2
""")

@test AoC.solve(day, 1, L[day]) == 138
@test AoC.solve(day, 2, L[day]) == 66
