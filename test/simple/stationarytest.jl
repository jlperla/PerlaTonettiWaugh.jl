# Generate default parameters. 
simple_params = @with_kw (γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5, ξ = 1.0, π_tilde = x -> 1)

# Test them
results_algebraic = stationary_algebraic_simple(simple_params());
@test results_algebraic.g ≈ 0.0211182826;
@test results_algebraic.ν ≈ 1.75369955156;
@test results_algebraic.v(0) ≈ 35.04962283;
@test results_algebraic.v(2) ≈ 165.31581267;
@test results_algebraic.v(5) ≈ 3312.7957099;

# Test for one particular grid. MATLAB solve_nonlinear_system = false. 
z1 = unique([linspace(0.0, 1.0, 500)' linspace(1.0, 5.0, 201)'])
results_num1 = stationary_numerical_simple(simple_params(), z1)
@test results_num1.g ≈ 0.020621554973703 atol = 1e-3 # Growth rate 

# Test for a third grid. 
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 11)' linspace(2.0, 5.0, 20)'])
results_num2 = stationary_numerical_simple(simple_params(), z)
@test results_num2.g ≈ 0.021404753523204 atol = 1e-2 # Growth rate 

# Test approximateness on the value function. 
ξ = simple_params().ξ
@test norm((results_num1.v.*exp.(ξ*z1))[1:end-7] - results_algebraic.v.(z1[1:end-7]), Inf) < 1 # Loosen this up. 

# Test for change zbar for grid and add points.
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 60)' linspace(2.0, 8.0, 40)'])
results = stationary_numerical_simple(simple_params(), z)
@test results.g ≈ 0.0211796240274 # Growth rate

# Test for change zbar for grid and add points.
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 60)' linspace(2.0, 10.0, 40)'])
results = stationary_numerical_simple(simple_params(), z)
@test results.g ≈ 0.02123967993879092 # Growth rate

# a baseline grid z
z = unique([linspace(0.0, 1.0, 300)' linspace(1.0, 2.0, 50)' linspace(2.0, 7.0, 50)'])
results = stationary_numerical_simple(simple_params(), z)
@test results.g ≈ 0.0211710310711 # Growth rate