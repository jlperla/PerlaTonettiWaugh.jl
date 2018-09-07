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
σ = baselineparams().σ

# Solve and compute residuals
@time solved = solve_dynamics(baselineparams(), settings(), d_0, d_T)

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