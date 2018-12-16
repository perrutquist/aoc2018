"""
    main(a::Symbol, v)

Assign `v` to a variable in the `Main` module, with name given by the symbol `a`,
and display it to stdout.

This can be useful for debugging, but should probably never be used in
production code.
"""
function main(a::Symbol, v, verbose)
    Main.eval( :( $a = $v ) )
    if verbose
        println(a, " =")
        display(v)
    end
end

"""
    @main x

Copy a variable (e.g. `x`) into the `Main` module, and display it to stdout.

This can be useful for debugging, but should probably never be used in
production code.
"""
macro main(a::Symbol, verbose::Bool=true)
    quote
        main($(QuoteNode(a)), $(esc(a)), $verbose)
    end
end
