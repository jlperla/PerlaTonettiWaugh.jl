using PerlaTonettiWaugh, Base.Test, Parameters, NamedTuples

# Generate default parameters. 

simple_numerical_params = @kw_nt(γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5)


z_min = 0
z_max = 5
M = 100
z = linspace(z_min, z_max, M);
temp = (rand(M) - 0.5)/100;
temp[1] = 0
temp[end] = 0
z = eye(M)*(z + temp);

# Test
results = stationary_numerical_simple(simple_numerical_params(), z)

@test results.g ≈ 0.0211182826;
@test results.ν ≈ 1.75369955156;
@test results.v[1] ≈ 35.04962283;
@test results.v[40] ≈ 165.31581267;
@test results.v[end] ≈ 3312.7957099;
