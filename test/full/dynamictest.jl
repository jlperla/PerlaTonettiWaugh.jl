# State grid. 
z_min = 0.0 
z_max = 5.0
M = 1000
z_grid = range(z_min, stop = z_max, length = M) # Since we only care about the grid. 

# Time
tstops_min_Δ_val = 1e-3 # minimum distance between tstops to be used for DE solvers.
# Define common objects. 
d_0 = 5
d_T = 2.3701
params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 
σ = params.σ
δ = params.δ

# Compute the stationary solution at t = 0 and t = T first
params_0 = merge(params, (d = d_0,)) # parameters to be used at t = 0
params_T = merge(params, (d = d_T,)) # parameters to be used at t = T

stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
Ω_0 = stationary_sol_0.Ω

stationary_sol_T = stationary_numerical(params_T, z_grid) # solution at t = T
v_T = stationary_sol_T.v_tilde
g_T = stationary_sol_T.g
z_hat_T = stationary_sol_T.z_hat
Ω_T = stationary_sol_T.Ω

# compute the resulting end time and function of Ω
T = (log(Ω_0) - log(Ω_T)) / δ

# T = sqrt(2*(log(Ω_0) - log(Ω_T)) / δ)
# Ω(t) = t < T ? Ω_0 * exp(-δ*T*t + δ*t^2/2) : Ω_T

T = 20.0
Ω(t) = Ω_T

settings = (z = z_grid, tstops = 0:1e-3:T, Δ_E = 1e-04)
