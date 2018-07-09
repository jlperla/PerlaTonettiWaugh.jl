#Generator for tuples with keyword arguments
macro kw_nt(args...)
    splits = map(args) do arg
        @match arg begin
            (a_ = b_) => (a, b)
            any_ => error("All arguments must be assignments")
        end
    end
    esc(:(
        (;$(map(splits) do pair
            Expr(:kw, pair[1], pair[2])
        end...),) -> 
        $NamedTuples.@NT($(map(splits) do pair
            Expr(:kw, pair[1], pair[1])
        end...))
    ))
end
