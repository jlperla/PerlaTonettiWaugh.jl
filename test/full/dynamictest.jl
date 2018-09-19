# User settings. 
  # Time and space grid. 
    z_min = 0.0 
    z_max = 5.0 
    M = 1000
    T = 20.0 
  # Experiment settings. 
    d_0 = 5.0
    d_T = 2.3701
    Δ_E = 1e-06
  # Overall parameters. 
    params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053)
  # Solver settings. 
    tstops = 0:1e-3:T # We don't currently use this anywhere. 

# Construct intermediate objects. 
  z = range(z_min, stop = z_max, length = M)
  settings = (z = z, tstops = tstops, Δ_E = Δ_E)
  params_0 = merge(params, (d = d_0,))
  params_T = merge(params, (d = d_T,))
  
# Get (numerical) stationary solution. 
  stationary_0 = stationary_numerical(params_0, z)
  stationary_T = stationary_numerical(params_T, z)

# Process those 
  Ω_0 = stationary_0.Ω
  Ω_T = stationary_T.Ω
  v_T = stationary_T.v_tilde
  g_T = stationary_T.g 
  z_hat_T = stationary_T.z_hat 

# Define more interim quantities. 
  Ω = t -> Ω_T # This is constant. 
  E = t -> (log(Ω(t + Δ_E)) - (log(t - Δ_E)))/(2*Δ_E) + params.δ # Log forward differences. 

@testset "Regression Tests" begin 
  # Run the solver. 
    sol = solve_dynamics(params_T, stationary_T, settings, T, Ω)
  # Spot-checks.
    @test sol.sol.t[5] ≈ 19.840000000000003
    @test sol.results[:λ_ii][end] ≈ 0.7813233366790822
    @test sol.sol.u[4][3] ≈ 1.1499093016834008
    @test sol.sol.prob.u0[1] ≈ 1.1868000000001158
  # Consistency checks. 
    @test norm(sol.results[:g] .- 0.0204198805) ≈ 0.0 atol = 1e-8 # Probably the most important of these checks. 
    @test norm(sol.results[:z_hat] .- 1.42571253) ≈ 0.0 atol = 1e-7 
    @test norm(sol.results[:Ω] .- 1.036549465138955) ≈ 0.0 atol = 1e-8
    @test all(sol.results[:E] .== 0.053)
    @test norm(sol.results[:S] .-  0.058475033) ≈ 0.0 atol = 1e-8 
    @test norm(sol.results[:L_tilde] .- 0.200430770) ≈ 0.0 atol = 1e-8

  # Run the solver for another case. 
    sol = solve_dynamics(params_0, stationary_0, settings, T, Ω)
  # Spot-checks. 
    @test sol.sol.t[5] == 19.852
    @test sol.results[:λ_ii][1] ≈ 0.9929472025880611
    @test sol.sol.u[4][3] ≈ 1.153063927522336
    @test sol.sol.prob.u0[1] ≈ 1.1868000000002454
  # Detailed checks. 
    @test sol.results[:g][1] ≈ 0.020019475192487802 # g check.
    @test sol.results[:g][end] ≈ 0.007963191154810903 # g check. 
    @test sol.results[:z_hat][1] ≈ 2.771561823423923 
    @test sol.results[:z_hat][(end-9)] ≈ 2.77021657056094 atol = 1e-8
end

@testset "Correctness Tests" begin # Here, we compare the DAE output to known correct values, such as MATLAB output or analytical results.
  # First case. 
  sol = solve_dynamics(params_T, stationary_T, settings, T, Ω)
  @test all([isapprox(x, 0.0, atol = 1e-9) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t. 

  # Second case. 
  sol = solve_dynamics(params_0, stationary_0, settings, T, Ω)
  # @show all([isapprox(x, 0.0, atol = 1e-6) for x in sol.results[:entry_residual]]) # Free-entry condition holds ∀ t.   
end 

@testset "Interpolation and entry_residuals Tests" begin 
  # Objects for interpolation. 
    Ω_nodes = 0:1e-1:T
    entry_residuals_nodes = Ω_nodes
    Ω_vec = map(t -> Ω(t), Ω_nodes)
  # First case. 
    @time sol = solve_dynamics(params_T, stationary_T, settings, T, Ω_vec, Ω_nodes)
  # Tests. 
    @test mean(sol.results[:entry_residual]) ≈ 0.0 atol = 1e-10
    residuals_interp = entry_residuals(params_T, stationary_T, settings, T, Ω_vec, Ω_nodes, entry_residuals_nodes).entry_residuals_interpolation
    @test mean(residuals_interp.(Ω_nodes)) ≈ 0.0 atol = 1e-9
end 

@testset "Regression Tests for f! at T" begin 
  # Instantiate the f! arguments
    du = zeros(M+2)
    u = [v_T..., g_T, z_hat_T]
    resid = zeros(M+2)
    t = T 
    p = []
  # Fill the residuals vector 
    f! = solve_dynamics(params_T, stationary_T, settings, T, Ω).f!
    f!(resid, du, u, p, t) # Kept the p from the old tests, since the solver complains without it. But it's a dummy. 
  # Tests
    # Accuracy
      @test mean(resid[1:M]) ≈ 0 atol = 1e-8
      @test mean(resid[M+1]) ≈ 0 atol = 1e-8
      @test mean(resid[M+2]) ≈ 0 atol = 1e-8
      @test all(abs.(resid) .<= 1e-6) # Test that we have small residuals across the board. 
    # Regression 
      @test resid[4] ≈ 0 atol = 1e-8
      @test resid[100] ≈ 0 atol = 1e-8 
end 