# Generate default parameters. 
simple_algebraic_params = @with_kw (γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5)

# Test them
results = stationary_algebraic_simple(simple_algebraic_params());
@test results.g ≈ 0.0211182826;
@test results.ν ≈ 1.75369955156;
@test results.v(0) ≈ 35.04962283;
@test results.v(2) ≈ 165.31581267;
@test results.v(5) ≈ 3312.7957099;
