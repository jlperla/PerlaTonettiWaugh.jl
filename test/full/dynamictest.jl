# User settings.
  # Time and space grid.
    z_min = 0.0
    z_max = 5.0
    P = 75
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

@testset "Regression Tests with constant Ω" begin
  # Run the solver.
    sol = solve_dynamics(params_T, stationary_T, settings, T, Ω, E)
  # Spot-checks.
    @test sol.sol.t[5] ≈ 19.840000000000003
    @test sol.results[:λ_ii][end] ≈ 0.8205696009690461
    @test sol.sol.u[4][3] ≈ 0.713179042868647
    @test sol.sol.prob.u0[1] ≈ 0.9329809580656181
  # Consistency checks.
    @test mean(sol.results[:g]) ≈ 0.006332516508877774 # Probably the most important of these checks.
    @test mean(sol.results[:z_hat]) ≈ 1.626173967797714
    @test mean(sol.results[:Ω]) ≈ 1.6928323073432825
    @test all(sol.results[:E] .== 0.053)
    @test mean(sol.results[:S]) ≈ -0.013749473107398985
    @test mean(sol.results[:L_tilde]) ≈ 0.1892988756097014
  # Full coverage of each column.
    @test sol.results[:t][4] == 1.5
    @test sol.results[:g][5] ≈ 0.006332516539866723
    @test sol.results[:z_hat][6] ≈ 1.6261739679041103
    @test sol.results[:Ω][7] ≈ 1.6928323073432825
    @test sol.results[:E][8] == 0.053
    @test sol.results[:v_1][9] ≈ 0.9329809601676073
    @test sol.results[:L_tilde][10] ≈ 0.18929887587415023
    @test sol.results[:λ_ii][11] ≈ 0.8205696004618322
    @test sol.results[:c][12] ≈ 0.9548509646788975
    @test sol.results[:S][4] ≈ -0.013749472948335131
    @test sol.results[:z_bar][3] ≈ 1.1778088573728906
    @test sol.results[:π_min][2] ≈ 0.14649422762537112
    @test sol.results[:entry_residual][12] ≈ 2.7166018323754315e-9 atol = 1e-5

  # TODO: Move this to the robustness section, and keep the original 1000-node grid (it doesn't work otherwise.)
  # # Run the solver for another case.
  #   sol = solve_dynamics(params_0, stationary_0, settings, T, Ω, E)
  # # Spot-checks.
  #   @test sol.sol.t[5] == 19.852
  #   @test_broken sol.results[:λ_ii][1] ≈ 0.9929472025880611
  #   @test_broken sol.sol.u[4][3] ≈ 1.153063927522336
  #   @test_broken sol.sol.prob.u0[1] ≈ 1.1868000000002454
  # # Detailed checks.
  #   @test_broken sol.results[:g][1] ≈ 0.020019475192487802 # g check.
  #   @test_broken sol.results[:z_hat][1] ≈ 2.771561823423923
  #   @test_broken sol.results[:z_hat][(end-9)] ≈ 2.770363670724641
  #   # Sub-pieces of L_tilde
  #   @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
  #   # Check if π_rat definitions in dynamics solution and SS coincide
  #   @test_broken sol.results[:π_rat][end] ≈ stationary_T.π_rat atol = 1e-3
end

@testset "Regression Tests with time-varying Ω" begin
    # Run the solver.
    T = sqrt(2*(log(Ω_0) - log(Ω_T)) / params.δ)
    Ω_t(t) = t < T ? Ω_0 * exp(-params.δ*T*t + params.δ*t*t/2) : Ω_T # Exponential Ω with time smoothing
    E_t(t) = (log(Ω_t(t + Δ_E)) - (log(Ω_t(t - Δ_E))))/(2*Δ_E) + params.δ # central differences.
    sol = solve_dynamics(params_T, stationary_T, settings, T, Ω_t, E_t)

    # Spot-checks.
    @test sol.sol.t[5] ≈ 3.532316848641289 atol = 1e-5
    @test sol.results[:λ_ii][end] ≈ 0.8205696009690461 atol = 1e-5
    @test sol.sol.u[4][3] ≈ 0.713186458707734 atol = 1e-5
    @test sol.sol.prob.u0[1] ≈ 0.9329809580656181atol = 1e-5
    # Consistency checks.
    @test mean(sol.results[:g]) ≈ 0.020144154416595162 atol = 1e-5 # Probably the most important of these checks.
    @test mean(sol.results[:z_hat]) ≈ 1.6118214712011016 atol = 1e-5
    @test mean(sol.results[:Ω]) ≈ 1.8255265725137522 atol = 1e-5
    @test mean(sol.results[:S]) ≈ 0.05706141328167729 atol = 1e-5
    @test mean(sol.results[:L_tilde]) ≈ 0.0926087309072632 atol = 1e-5
    # Full coverage of each column.
    @test sol.results[:t][4] ≈ 0.8318036450026267 atol = 1e-5
    @test sol.results[:g][5] ≈ 0.03529457859081019 atol = 1e-5
    @test sol.results[:z_hat][6] ≈ 1.5963201866789414 atol = 1e-5
    @test sol.results[:Ω][7] ≈ 1.7902142064165565 atol = 1e-5
    @test sol.results[:v_1][9] ≈ 0.98515301205966 atol = 1e-5
    @test sol.results[:L_tilde][10] ≈ 0.160812673571405 atol = 1e-5
    @test sol.results[:λ_ii][11] ≈ 0.8197747097913394 atol = 1e-5
    @test sol.results[:c][12] ≈ 0.9667053570460247 atol = 1e-5
    @test sol.results[:S][4] ≈ 0.14915482412907646 atol = 1e-5
    @test sol.results[:z_bar][3] ≈ 1.2083652214290672 atol = 1e-5
    @test sol.results[:π_min][2] ≈ 0.2042841017741823 atol = 1e-5
    @test sol.results[:entry_residual][12] ≈ 0.013589358773159477 atol = 1e-5
    # Sub-pieces of L_tilde
    @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
    # TODO: is this still valid with a coarse grid? 
    # Check if π_rat definitions in dynamics solution and SS coincide
    # @test_broken sol.results[:π_rat][end] ≈ stationary_T.π_rat atol = 1e-3
end

@testset "Correctness Tests" begin # Here, we compare the DAE output to known correct values, such as MATLAB output or analytical results.
  # First case.
  sol = solve_dynamics(params_T, stationary_T, settings, T, Ω, E)
  @test all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.

  # TODO: Move this to robustness, with the original fine grid.
  # sol = solve_dynamics(params_0, stationary_0, settings, T, Ω, E)
  # @test all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.
end

@testset "Regression Tests for f! at T" begin
  # Instantiate the f! arguments
    du = zeros(P+2)
    u = [v_T..., g_T, z_hat_T]
    resid = zeros(P+2)
    t = T
    p = solve_dynamics(params_T, stationary_T, settings, T, Ω, E).p
  # Fill the residuals vector
    f!(resid, du, u, p, t)
  # Tests
    # Accuracy
      @test mean(resid[1:P]) ≈ 0 atol = 1e-6
      @test mean(resid[P+1]) ≈ 0 atol = 1e-6
      @test mean(resid[P+2]) ≈ 0 atol = 1e-6
      @test all(abs.(resid) .<= 1e-6) # Test that we have small residuals across the board.
    # Regression
      @test resid[4] ≈ 0 atol = 1e-6
      @test resid[end] ≈ 0 atol = 1e-6
    # Benchmarking results
      println("Kernel Benchmarking")
      @show @btime f!($resid, $du, $u, $p, $t)
end
