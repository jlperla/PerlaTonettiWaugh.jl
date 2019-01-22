# Generate default parameters.
simple_params = @with_kw (μ = 0.0048, υ = 0.02, θ = 2.1, r = 0.05, ζ = 14.5, ξ = 1.0, π = x -> 1)

# Test them
results_algebraic = stationary_algebraic_simple(simple_params());
@test results_algebraic.g ≈ 0.020727753490265812; # from updating the eq in stationary
@test results_algebraic.ν ≈ 1.797254207263066;
@test results_algebraic.v(0) ≈ 34.58676260411018;
@test results_algebraic.v(2)*exp(2) ≈ 164.54095232228548; # correct for rescaling
@test results_algebraic.v(5)*exp(5) ≈ 3298.0717492785516;

# Test for one particular grid.
z1 = unique([range(0.0, stop = 1.0, length = 500)' range(1.0, stop = 5.0, length = 201)'])
results_num1 = stationary_numerical_simple(simple_params(), z1)
@test results_num1.g ≈ 0.02010032283241259 # Invariance
@test abs(results_num1.g - results_algebraic.g) < 1e-3 # Consistency.

# Test for a second grid.
z = unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 11)' range(2.0, stop = 5.0, length = 20)'])
results_num2 = stationary_numerical_simple(simple_params(), z)
@test abs(results_num2.g - results_algebraic.g) < 1e-3 # Growth rate

# Test approximateness on the value function.
ξ = simple_params().ξ
@test norm((results_num2.v.*exp.(ξ*z))[1:end-7] - results_algebraic.v.(z[1:end-7]) .* exp.(ξ*z[1:end-7]), Inf) < 1 # Tighten this up.

# Test for change zbar for grid and add points.
z = unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 60)' range(2.0, stop = 8.0, length = 40)'])
results = stationary_numerical_simple(simple_params(), z)
@test abs(results.g - results_algebraic.g) < 1e-4  # Growth rate

# Test for change zbar for grid and add points.
z = unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 60)' range(2.0, stop = 10.0, length = 40)'])
results = stationary_numerical_simple(simple_params(), z)
@test abs(results.g - results_num1.g) < 1e-3 # Growth rate

# a baseline grid z
z = unique([range(0.0, stop = 1.0, length = 300)' range(1.0, stop = 2.0, length = 50)' range(2.0, stop = 7.0, length = 50)'])
results = stationary_numerical_simple(simple_params(), z)
@test abs(results.g - results_algebraic.g) < 1e-4  # Growth rate
