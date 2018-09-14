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

# Solve and compute residuals
@time solved = solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)

@test mean(mean(solved.residuals[:,1:M], dims = 1)) ≈ 0 atol = 1e-03 # mean residuals for system of ODEs
@test mean(mean(solved.residuals[:,(M+1)])) ≈ 0 atol = 1e-03 # mean residuals for value matching condition
@test mean(mean(solved.residuals[:,(M+2)])) ≈ 0 atol = 1e-03 # mean residuals for export threshold condition

v_hat_t0 = map(z -> exp((σ-1)*z), z_grid) .* solved.v[1]
@test any((v -> v < 0).(diff(v_hat_t0))) == false # after reparametrization v_hat should be increasing

# # plot v_hat, v, g, and z_hat in `test/full/ptw_plots/`
# using Plots
# v0 = map(v -> v[1], solved.v)
# plot(z_grid, v_hat_t0, label = "v_hat at t = 0", lw = 3)
# savefig("test/full/ptw_plots/v_hat_t0.png")
# plot(z_grid, solved.v[1], label = "v at t = 0", lw = 3)
# savefig("test/full/ptw_plots/v_t0.png")
# plot(z_grid, solved.v[end], label = "v at t = T", lw = 3)
# savefig("test/full/ptw_plots/v_tT.png")
# plot(solved.t, v0, label = "v0", lw = 3)
# savefig("test/full/ptw_plots/v0.png")
# plot(solved.t, solved.g, label = "g", lw = 3)
# savefig("test/full/ptw_plots/g.png")
# plot(solved.t, solved.z_hat, label = "z_hat", lw = 3)
# savefig("test/full/ptw_plots/z_hat.png")
# # ... and S, L_tilde, z_bar, π_min, π_tilde_z0, λ_ii, c, E
# plot(solved.t, solved.S, label = "S", lw = 3)
# savefig("test/full/ptw_plots/S.png")
# plot(solved.t, solved.L_tilde, label = "L_tilde", lw = 3)
# savefig("test/full/ptw_plots/L_tilde.png")
# plot(solved.t, solved.z_bar, label = "z_bar", lw = 3)
# savefig("test/full/ptw_plots/z_bar.png")
# plot(solved.t, solved.π_min, label = "pi_min", lw = 3)
# savefig("test/full/ptw_plots/π_min.png")
# plot(z_grid, solved.π_tilde[1], label = "pi_tilde at t = 0", lw = 3)
# savefig("test/full/ptw_plots/π_tilde_t0.png")
# plot(solved.t, solved.λ_ii, label = "lambda_ii", lw = 3)
# savefig("test/full/ptw_plots/λ_ii.png")
# plot(solved.t, solved.c, label = "c", lw = 3)
# savefig("test/full/ptw_plots/c.png")
# plot(solved.t, solved.E, label = "E", lw = 3)
# savefig("test/full/ptw_plots/E.png")
# # ... and entry_residual
# plot(solved.t, solved.entry_residual, label = "entry_residual", lw = 3)
# savefig("test/full/ptw_plots/entry_residual.png")