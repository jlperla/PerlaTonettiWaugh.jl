include("matlabobjects.jl")

# Generate default parameters. 
simple_numerical_params = @with_kw (γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5, ξ=1)

# Test for one particular grid. MATLAB solve_nonlinear_system = false. 
z = unique([linspace(0.0, 1.0, 500)' linspace(1.0, 5.0, 201)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test results.g ≈ 0.0204899363906 # Growth rate 

# Test for a third grid. 
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 11)' linspace(2.0, 5.0, 20)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test results.g ≈ 0.02058576255994 # Growth rate 

# Test for change zbar for grid and add points. 
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 60)' linspace(2.0, 8.0, 40)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test results.g ≈ 0.0211796240274 # Growth rate 
