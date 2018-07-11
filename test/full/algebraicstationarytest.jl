@suppress using PerlaTonettiWaugh, Base.Test, NLsolve, BenchmarkTools
include("../../src/full/fullparams.jl")

# Run a basic test. (neglecting nlsolve stability)
MATLAB_eq = [0.0260, 3.7504, 1.0615] # Equilibrium for default fullparams() from MATLAB.
init_x = copy(MATLAB_eq)
@btime result = nlsolve(f!, init_x)
@btime result = nlsolve(f!, [0.01, 4.3, 2.0])

result = nlsolve(f!, init_x)
@test round.(result.zero, 4) ≈ MATLAB_eq # Test if we're nearby 
result = nlsolve(f!, [0.01, 4.3, 2.0])
@test round.(result.zero, 4) ≈ MATLAB_eq # Test success for a different starting point. 

@test_throws AssertionError nlsolve(f!, [0.01, 2.2, 2.0]) # Test for <0 value failure (low z_hat guess.)
@test_throws AssertionError nlsolve(f!, [0.01, 2.2, 1.0]) # Another broken test.  
@test_throws AssertionError nlsolve(f!, [0.01, 2.2, 2.0]; method = :anderson) # Test failure for Anderson.

# Different case (custom params test)
params2 = @kw_nt(ρ = 0.5, σ = 3.0, N = 10.0, θ = 3.22, γ = 0, d = 5.49, κ = 0.06, ζ = 1.9, η = 0, Theta = 1.0, χ = 1/3, υ = 0.01, μ = 0, δ = 0.01)
result = nlsolve((G, x) -> f!(G, x; params=params2), [0.12, 1.2, 0.04])
@test round.(result.zero, 5) ≈ [0.24444, 1.33895, 0.12194]

# Benchmarking
init_x = [0.01, 6.0, 0.1]
@btime nlsolve(f!, init_x)