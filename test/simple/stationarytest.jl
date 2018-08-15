# Generate default parameters. 
simple_params = @with_kw (μ = 0.0048, υ = 0.02, θ = 2.1, r = 0.05, ζ = 14.5, ξ = 1.0, π_tilde = x -> 1)

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


#= 
baselineparams = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 
=#

baselineparams = simple_params(θ = 5.1269, υ = 0.0593, μ = 0.0, ζ = 1) # μ should probably be 0.0
results_baseline_algebraic = stationary_algebraic_simple(baselineparams)
results_baseline_numerical = stationary_numerical_simple(baselineparams, z)
@show diffsols(results_baseline_numerical, results_baseline_algebraic)