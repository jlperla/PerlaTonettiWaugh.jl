using Base.Test, PerlaTonettiWaugh 

# Define parameters and get results from each function. 
simple_algebraic_params = @with_kw(γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5)
algebraic = stationary_algebraic_simple(simple_algebraic_params());
# numerical = stationary_numerical_simple(simple_algebraic_params()); # Or whatever the API is. 

# Test for equality 
#= 
    @test algebraic.g ≈ numerical.g
    @test algebraic.ν ≈ numerical.ν
    @test algebraic.v(0) ≈ numerical.v(0)
    @test algebraic.v(2) ≈ numerical.v(2)
    @test algebraic.v(5) ≈ numerical.v(5)
    # Something for the other objects? 
=# 