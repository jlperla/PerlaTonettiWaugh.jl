# NOTE: This file is mainly for records-keeping purposes, and the interface may have evolved.
# User settings.
  # Time and space grid.
  z_min = 0.0
  z_max = 5.0
  P = 1000
  T = 20.0
# Experiment settings.
  d_0 = 5.0
  d_T = 2.3701
  Δ_E = 1e-06
# Overall parameters.
  params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053)
# Solver settings.
  tstops = nothing

# Construct intermediate objects.
z_ex = range(z_min, stop = z_max, length = P+2) 
z = z_ex[2:end-1]
settings = (z = z, z_ex = z_ex, tstops = tstops, Δ_E = Δ_E, T_U_bar = 1.0)
params_0 = merge(params, (d = d_0,))
params_T = merge(params, (d = d_T,))

# Get (numerical) stationary solution.
stationary_0 = stationary_numerical(params_0, z_ex)
stationary_T = stationary_numerical(params_T, z_ex)

# Process those
Ω_0 = stationary_0.Ω
Ω_T = stationary_T.Ω
v_T = stationary_T.v_tilde
g_T = stationary_T.g
z_hat_T = stationary_T.z_hat

# Define more interim quantities.
Ω(t) = Ω_T # This is constant.
E(t) = (log(Ω(t + Δ_E)) - (log(Ω(t - Δ_E))))/(2*Δ_E) + params.δ # central differences.

@testset "Constant Omega" begin 
    sol = solve_dynamics(params_0, stationary_0, settings, T, Ω, E)
  # Spot-checks.
    @test sol.sol.t[5] == 19.852
    @test sol.results[:λ_ii][1] ≈ 0.993204827349257
    @test sol.sol.u[4][3] ≈ 1.1338849871884573 
    @test sol.sol.prob.u0[1] ≈ 1.1675290237815987
  # Detailed checks.
    @test sol.results[:g][1] ≈ 0.01739899362603311 # g check.
    @test sol.results[:z_hat][1] ≈  2.827474674128036
    @test sol.results[:z_hat][(end-9)] ≈ 2.8259905922985666
    # Sub-pieces of L_tilde
    @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
end
