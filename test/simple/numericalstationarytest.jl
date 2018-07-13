using PerlaTonettiWaugh, Base.Test, Parameters, NamedTuples
include("matlabobjects.jl")
# Generate default parameters. 

simple_numerical_params = @kw_nt(γ = 0.005, σ = 0.02, α = 2.1, r = 0.05, ζ = 14.5)

# Test for one particular grid. MATLAB solve_nonlinear_system = false. 
z = unique([linspace(0.0, 1.0, 500)' linspace(1.0, 5.0, 201)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test results.g ≈ 0.020621554973703 # Growth rate 
@test_broken results.ν ≈ 1.753699551569504 # Nu 
@test results.v ≈ val1 # Value Function

# Test for another grid. 
z = unique([linspace(0.0, 1.0, 50)' linspace(1.0, 2.0, 20)' linspace(1.0, 5.0, 500)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test_broken results.g ≈ 0.021020961277543 # Growth rate
@test_broken results.ν ≈ 1.753699551569504 # Nu
@test_broken results.v ≈ val2 # Value function

# Test for a third grid. 
z = unique([linspace(0.0, 1.0, 1000)' linspace(1.0, 2.0, 11)' linspace(2.0, 5.0, 20)'])
results = stationary_numerical_simple(simple_numerical_params(), z)
@test results.g ≈ 0.021404753523204 # Growth rate 
@test results.ν ≈ 1.7075596085501772 # Nu 
@test results.v ≈ val3 # Value function
