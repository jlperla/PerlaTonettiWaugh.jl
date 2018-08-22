using PerlaTonettiWaugh, Base.Test
using Distributions, Sundials, BenchmarkTools, QuantEcon, Interpolations, Parameters, NamedTuples, NLsolve, ContinuousTransformations, DifferentialEquations

# State grid. 
z_min = 0.0 
z_max = 5.0
M = 100
z_grid = linspace(z_min,z_max,M) # Since we only care about the grid. 

# Time
T_val = 100.0

# Define common objects. 
baselineparams = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 
settings = @with_kw (z = z_grid, T = T_val)

solved = solve_dynamic_full(baselineparams(), settings())