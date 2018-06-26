# Unit testing for the algebraic stationary solution, simple model. 

# Arnav Sood: Jun. 26, 2018

using Base.Test, NamedTuples, Parameters, MacroTools
include("..\\..\\src\\simple\\algebraicstationary.jl"); 

#Generator for tuples with keyw
macro with_kw(args...)
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

# Generate default parameters. 
params = @with_kw(γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5)

# Test them
results = stationary_algebraic(params());
@test results.g ≈ 0.0211182826;
@test results.υ ≈ 1.75369955156;
@test results.v(0) ≈ 35.04962283;
@test results.v(2) ≈ 165.31581267;
@test results.v(5) ≈ 3312.7957099;
