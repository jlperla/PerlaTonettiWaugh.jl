# User settings.
  # Time and space grid.
    z_min = 0.0
    z_max = 5.0
    P = 75
    T_default = 20.0
  # Experiment settings.
    d_0 = 5.0
    d_T = 2.3701
    Δ_E = 1e-06
  # Overall parameters.
    params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053)
  # Solver settings.

# Construct intermediate objects.
  z_ex = range(z_min, stop = z_max, length = P+2)
  z = z_ex[2:end-1]
  settings = settings_defaults(z_ex = z_ex, T = T_default)
  params_0 = merge(params, (d = d_0,))
  params_T = merge(params, (d = d_T,))

# Get (numerical) stationary solution.
  stationary_0 = stationary_numerical(params_0, settings_defaults(z_ex = z_ex))
  stationary_T = stationary_numerical(params_T, settings_defaults(z_ex = z_ex))

# Process those
  Ω_0 = stationary_0.Ω
  Ω_T = stationary_T.Ω
  v_T = stationary_T.v_tilde
  g_T = stationary_T.g
  z_hat_T = stationary_T.z_hat

# Define more interim quantities.
  Ω(t) = Ω_T # This is constant.

@testset "Regression Tests with constant Ω" begin
  # Run the solver.
    sol = solve_dynamics(params_T, stationary_T, settings, Ω)
    sol = merge(sol, (results = prepare_results(sol, stationary_T, stationary_0),))
    filter!(row -> row.t >= 0, sol.results) # get rid of pre-shock values 

    # Spot-checks.
    @test sol.sol.t[5] ≈ 19.996000000000002 atol = 1e-5
    @test sol.results[:λ_ii][end] ≈ 0.8205696009690461 atol = 1e-5
    @test sol.sol.u[4][3] ≈ 0.713179042868647 atol = 1e-5
    @test sol.sol.prob.u0[1] ≈ 0.9329809580656181 atol = 1e-5
  # Consistency checks.
    @test mean(sol.results[:g]) ≈ 0.006332516508877774  atol = 1e-5 # Probably the most important of these checks.
    @test mean(sol.results[:z_hat]) ≈ 1.626173967797714 atol = 1e-5
    @test mean(sol.results[:Ω]) ≈ 1.6928323073432825 atol = 1e-5
    # @test all(sol.results[:E] .== 0.053)
    @test mean(sol.results[:S]) ≈ -0.013749473107398985 atol = 1e-5
    @test mean(sol.results[:L_tilde]) ≈ 0.1892988756097014 atol = 1e-5
  # Full coverage of each column.
    @test sol.results[:t][4] ≈ 0.75 atol = 1e-5
    @test sol.results[:g][5] ≈ 0.006332516539866723 atol = 1e-5
    @test sol.results[:z_hat][6] ≈ 1.6261739679041103 atol = 1e-5
    @test sol.results[:Ω][7] ≈ 1.6928323073432825 atol = 1e-5
    @test sol.results[:E][8] ≈ 0.05300001047668845 atol = 1e-5
    @test sol.results[:v_1][9] ≈ 0.9329809601676073 atol = 1e-5
    @test sol.results[:L_tilde][10] ≈ 0.18929887587415023 atol = 1e-5
    @test sol.results[:λ_ii][11] ≈ 0.8205696004618322 atol = 1e-5
    @test sol.results[:c][12] ≈ 1.3801007368941405 atol = 1e-5
    @test sol.results[:S][4] ≈ -0.013749472948335131 atol = 1e-5
    @test sol.results[:z_bar][3] ≈ 1.702354537311347 atol = 1e-5
    @test sol.results[:π_min][2] ≈ 0.04423572434773227 atol = 1e-5
    @test sol.results[:entry_residual][12] ≈ 2.7166018323754315e-9 atol = 1e-5
    @test sol.results[:r][19] ≈ 0.0793325164566891 atol = 1e-5
end

@testset "Regression Tests with time-varying Ω" begin
    # Run the solver.
    T = sqrt(2*(log(Ω_0) - log(Ω_T)) / params.δ)
    Ω_t(t) = t < T ? Ω_0 * exp(-params.δ*T*t + params.δ*t*t/2) : Ω_T # Exponential Ω with time smoothing
    settings = settings_defaults(z_ex = z_ex, T = T)
    sol = solve_dynamics(params_T, stationary_T, settings, Ω_t)
    sol = merge(sol, (results = prepare_results(sol, stationary_T, stationary_0),))
    filter!(row -> row.t >= 0, sol.results)

    # Spot-checks.
    @test sol.sol.t[5] ≈ 3.77533130102673 atol = 1e-5
    @test sol.results[:λ_ii][end] ≈ 0.8205696009690461 atol = 1e-5
    @test sol.sol.u[4][3] ≈ 0.713186458707734 atol = 1e-5
    @test sol.sol.prob.u0[1] ≈ 0.9329809580656181 atol = 1e-5
    # Consistency checks.
    @test mean(sol.results[:g]) ≈ 0.005968185219011899 atol = 1e-5 # Probably the most important of these checks.
    @test mean(sol.results[:z_hat]) ≈ 1.6551587501780332 atol = 1e-5
    @test mean(sol.results[:Ω]) ≈ 1.8166904089164235 atol = 1e-5
    @test mean(sol.results[:S]) ≈ -0.01561736319741233 atol = 1e-5
    @test mean(sol.results[:L_tilde]) ≈ 0.18682781261445597 atol = 1e-5
    # Full coverage of each column.
    @test sol.results[:t][4] ≈ 0.75 atol = 1e-5
    @test sol.results[:g][5] ≈ 0.005164669443013719 atol = 1e-5
    @test sol.results[:z_hat][6] ≈ 1.7014106198547416 atol = 1e-5
    @test sol.results[:Ω][7] ≈ 1.9418121379049862 atol = 1e-5
    @test sol.results[:v_1][9] ≈ 0.9329809578947367 atol = 1e-5
    @test sol.results[:L_tilde][10] ≈ 0.18787808053459026 atol = 1e-5
    @test sol.results[:λ_ii][11] ≈ 0.8237900176167677 atol = 1e-5
    @test sol.results[:c][12] ≈ 1.3916274048069284 atol = 1e-5
    @test sol.results[:S][4] ≈ -0.020681135472245087 atol = 1e-5
    @test sol.results[:z_bar][3] ≈ 1.844294549955518 atol = 1e-5
    @test sol.results[:π_min][2] ≈ 0.03330744238041946 atol = 1e-5
    @test sol.results[:entry_residual][12] ≈ 0.0 atol = 1e-5
    @test sol.results[:r][19] ≈ 0.07705459878582416 atol = 1e-5
    # Sub-pieces of L_tilde
    @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
end

@testset "Correctness Tests" begin # Here, we compare the DAE output to known correct values, such as MATLAB output or analytical results.
  # First case.
  sol = solve_dynamics(params_T, stationary_T, settings, Ω)
  sol = merge(sol, (results = prepare_results(sol, stationary_T, stationary_0),))
  filter!(row -> row.t >= 0, sol.results)
  @test all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.
end

@testset "Regression Tests for f! at T" begin
  # Instantiate the f! arguments
    du = zeros(P+3)
    u = [v_T..., g_T, z_hat_T, params.δ]
    resid = zeros(P+3)
    t = T_default
    p = solve_dynamics(params_T, stationary_T, settings, Ω).dae_parameters
  # Fill the residuals vector
    f!(resid, du, u, p, t)
  # Tests
    # Accuracy
      @test mean(resid[1:P]) ≈ 0 atol = 1e-6
      @test mean(resid[P+1]) ≈ 0 atol = 1e-6
      @test mean(resid[P+2]) ≈ 0 atol = 1e-6
      @test mean(resid[P+3]) ≈ 0 atol = 1e-6
      @test all(abs.(resid) .<= 1e-6) # Test that we have small residuals across the board.
    # Regression
      @test resid[4] ≈ 0 atol = 1e-6
      @test resid[end] ≈ 0 atol = 1e-6
    # Benchmarking results
      println("Kernel Benchmarking")
      @show @btime f!($resid, $du, $u, $p, $t)
end
