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
@testset "Regression Tests with constant Ω" begin
  # Run the solver.
    sol = solve_dynamics(params_T, stationary_T, settings, T, Ω, E)
  # Spot-checks.
    @test sol.sol.t[5] ≈ 19.840000000000003
    @test_broken sol.results[:λ_ii][end] ≈ 0.7813233366790822
    @test_broken sol.sol.u[4][3] ≈ 1.1499093016834008
    @test sol.sol.prob.u0[1] ≈ 1.1868000000001158
  # Consistency checks.
    @test_broken mean(sol.results[:g]) ≈ 0.020419880565858414 # Probably the most important of these checks.
    @test_broken mean(sol.results[:z_hat]) ≈ 1.425712535030436
    @test_broken mean(sol.results[:Ω]) ≈ 1.036549465116929
    @test all(sol.results[:E] .== 0.053)
    @test_broken mean(sol.results[:S]) ≈ 0.05847503367633505
    @test_broken mean(sol.results[:L_tilde]) ≈ 0.200430770186847
  # Full coverage of each column.
    @test sol.results[:t][4] == 1.5
    @test_broken sol.results[:g][5] ≈ 0.020419880572700753
    @test_broken sol.results[:z_hat][6] ≈ 1.4257125350016342
    @test_broken sol.results[:Ω][7] ≈ 1.036549465116929
    @test sol.results[:E][8] == 0.053
    @test sol.results[:v_1][9] ≈ 1.1867999999850465
    @test_broken sol.results[:L_tilde][10] ≈ 0.20043077010213567
    @test_broken sol.results[:λ_ii][11] ≈ 0.7813233366545765
    @test_broken sol.results[:c][12] ≈ 1.1882939932974685
    @test_broken sol.results[:S][4] ≈ 0.05847503375286956
    @test_broken sol.results[:z_bar][3] ≈ 1.4861677378016562
    @test_broken sol.results[:π_min][2] ≈ 0.16549999774577948
    @test sol.results[:entry_residual][12] ≈ -1.3211653993039363e-13 atol = 1e-5

  # Run the solver for another case.
    sol = solve_dynamics(params_0, stationary_0, settings, T, Ω, E)
  # Spot-checks.
    @test sol.sol.t[5] == 19.852
    @test_broken sol.results[:λ_ii][1] ≈ 0.9929472025880611
    @test_broken sol.sol.u[4][3] ≈ 1.153063927522336
    @test sol.sol.prob.u0[1] ≈ 1.1868000000002454
  # Detailed checks.
    @test_broken sol.results[:g][1] ≈ 0.020019475192487802 # g check.
    @test_broken sol.results[:z_hat][1] ≈ 2.771561823423923
    @test_broken sol.results[:z_hat][(end-9)] ≈ 2.770363670724641
    # Sub-pieces of L_tilde
    @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
    # Check if π_rat definitions in dynamics solution and SS coincide
    @test_broken sol.results[:π_rat][end] ≈ stationary_T.π_rat atol = 1e-3
end

@testset "Regression Tests with time-varying Ω" begin
    # Run the solver.
    T = sqrt(2*(log(Ω_0) - log(Ω_T)) / params.δ)
    Ω_t(t) = t < T ? Ω_0 * exp(-params.δ*T*t + params.δ*t*t/2) : Ω_T # Exponential Ω with time smoothing
    E_t(t) = (log(Ω_t(t + Δ_E)) - (log(Ω_t(t - Δ_E))))/(2*Δ_E) + params.δ # central differences.
    sol = solve_dynamics(params_T, stationary_T, settings, T, Ω_t, E_t)

    # Spot-checks.
    @test_broken sol.sol.t[5] ≈ 4.0820840627115444 atol = 1e-5
    @test_broken sol.results[:λ_ii][end] ≈ 0.7813233366790822 atol = 1e-5
    @test_broken sol.sol.u[4][3] ≈ 1.1499297704103377 atol = 1e-5
    @test sol.sol.prob.u0[1] ≈ 1.1868000000001158 atol = 1e-5
    # Consistency checks.
    @test_broken mean(sol.results[:g]) ≈ 0.032776838190899334 atol = 1e-5 # Probably the most important of these checks.
    @test_broken mean(sol.results[:z_hat]) ≈ 1.4240459625635837 atol = 1e-5
    @test_broken mean(sol.results[:Ω]) ≈ 1.1546007592192749 atol = 1e-5
    @test_broken mean(sol.results[:S]) ≈ 0.12182791761996183 atol = 1e-5
    @test_broken mean(sol.results[:L_tilde]) ≈ 0.1068149526037354 atol = 1e-5
    # Full coverage of each column.
    @test_broken sol.results[:t][4] ≈ 1.2180412122607018 atol = 1e-5
    @test_broken sol.results[:g][5] ≈ 0.04344370969889273 atol = 1e-5
    @test_broken sol.results[:z_hat][6] ≈ 1.4138759332401472 atol = 1e-5
    @test_broken sol.results[:Ω][7] ≈ 1.0852902581449317 atol = 1e-5
    @test_broken sol.results[:v_1][9] ≈ 1.2162498772031762 atol = 1e-5
    @test_broken sol.results[:L_tilde][10] ≈ 0.18748845694858562 atol = 1e-5
    @test_broken sol.results[:λ_ii][11] ≈ 0.7810856222518524 atol = 1e-5
    @test_broken sol.results[:c][12] ≈ 1.1899458996247958 atol = 1e-5
    @test_broken sol.results[:S][4] ≈ 0.18806483288417494 atol = 1e-5
    @test_broken sol.results[:z_bar][3] ≈ 1.6350019300336807 atol = 1e-5
    @test_broken sol.results[:π_min][2] ≈ 0.21425367016810465 atol = 1e-5
    @test_broken sol.results[:entry_residual][12] ≈ 0.0010887433781643363 atol = 1e-5
    # Sub-pieces of L_tilde
    @test sol.results[:L_tilde_a] + sol.results[:L_tilde_x] + sol.results[:L_tilde_E] ≈ sol.results[:L_tilde]
    # Check if π_rat definitions in dynamics solution and SS coincide
    @test_broken sol.results[:π_rat][end] ≈ stationary_T.π_rat atol = 1e-3
end

@testset "Correctness Tests" begin # Here, we compare the DAE output to known correct values, such as MATLAB output or analytical results.
  # First case.
  sol = solve_dynamics(params_T, stationary_T, settings, T, Ω, E)
  @test all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.

  # Second case.
  sol = solve_dynamics(params_0, stationary_0, settings, T, Ω, E)
  # @show all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.
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
      @test resid[100] ≈ 0 atol = 1e-6
    # Benchmarking results
      println("Kernel Benchmarking")
      @show @btime f!($resid, $du, $u, $p, $t)
end
