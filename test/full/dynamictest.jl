# State grid. 
z_min = 0.0 
z_max = 5.0
M = 1000
z_grid = range(z_min, stop = z_max, length = M) # Since we only care about the grid. 

# Time
T_val = 100.0

# Define common objects. 
d_0 = 5
d_T = 2.3701
baselineparams = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 
settings = @with_kw (z = z_grid, T = T_val)

# Solve and compute residuals
@time solved = solve_dynamics(baselineparams(), settings(), d_0, d_T)

@test mean(mean(solved.residuals[:,1:M], dims = 1)) ≈ 0 atol = 1e-03 # mean residuals for system of ODEs
@test mean(mean(solved.residuals[:,(M+1)])) ≈ 0 atol = 1e-03 # mean residuals for value matching condition
@test mean(mean(solved.residuals[:,(M+2)])) ≈ 0 atol = 1e-03 # mean residuals for export threshold condition
