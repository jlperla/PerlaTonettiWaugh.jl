@testset "Simple Model" begin
# UTILITIES
# Diffusion operators
function testdiffusion()
    ξ = 1.0
    z = 0.:0.01:5.
    rescaled_diffusionoperators(z, ξ), rescaled_diffusionoperators(collect(z), ξ) # test both the AbstractRange and StepRangeLen methods
end

@test @inferred testdiffusion() == testdiffusion()

# ω_weights
function test_ω()
    z = 0.:0.01:5.
    α = 2.1
    ξ = 1.
    ω_weights(z, α, ξ), ω_weights(collect(z), α, ξ), ω_weights(z, α, 1), ω_weights(collect(z), α, 1), ω_weights(z, 1.0, 2.0) # Different type signatures/values for input
end

@test @inferred test_ω() == test_ω()

# DYNAMIC OBJECTS
# simpleODE
function test_simpleODE()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10., g = x -> 0.02)
    simpleODE(params, settings)
end

@btime test_simpleODE()

# simpleDAE
function test_simpleDAE()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10., g = x -> 0.02)
    simpleDAE(params, settings)
end

@btime test_simpleDAE()

# RESIDUALS CALCULATION
function test_residuals()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10, g = x -> 0.02, ode_solve_algorithm = CVODE_BDF(), t_grid = range(0., stop = 10, length = length(0.:0.1:3.)))
    obj1 = calculate_residuals(params, settings)

    # Vary the parameters a bit
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1. + 0.1*x, x = n -> 1, ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10, g = x -> 0.02, ode_solve_algorithm = CVODE_BDF(), t_grid = range(0., stop = 10, length = length(0.:0.1:3.)))
    obj2 = calculate_residuals(params, settings) # calls second method

    return obj1, obj2
end

@test @inferred test_residuals() == test_residuals()
end

@testset "Full Model" begin
# SETUP OBJECTS.
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
  tstops = nothing
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
# Arguments.
  du = zeros(M+2)
  u = [v_T..., g_T, z_hat_T]
  resid = zeros(M+2)
  t = T
  p = solve_dynamics(params_T, stationary_T, settings, T, Ω, E).p
# Test inference.
  # Out-of-place function.
  function f(du,u,p,t)
    # Setup (unpack arguments, reset residual, grab E and Ω evaluations, etc.)
      @unpack ζ, Ω, E, static_equilibrium, T, results, ρ, δ, σ, μ, υ, L_1, L_2, ω, κ, d = p
      residual .= 0
      M = length(residual) - 2
      # v = u[1:M] (this line is commented to cut down on memory allocations)
      g = u[M+1]
      z_hat = u[M+2]
      x = ζ
      Ω_t = Ω(t)
      E_t = E(t)
    # Get static equilibrium values
      @unpack S_t, L_tilde_t, z_bar, π_min, π_tilde = static_equilibrium(u[1], g, z_hat, E_t, Ω_t)
    # Grab the L_tilde derivative.
      L_tilde_log_derivative = 0.0 # Default to Float literal.
      if (t < T)
        t_forward = results[:t][end]
        L_tilde_forward = results[:L_tilde][end]
        L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde_t))/(t_forward - t) # Forward differences.
      end
    #=  Reset the residuals to slack in the DAE conditions.
        Note that A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2 and we're decomposing this.
    =#
      residual = fill(0.0, M+2)
      residual[1:M] = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:M] # system of ODEs (eq:28)
      residual[1:M] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:M]
      residual[1:M] .-= (υ^2/2)*L_2*u[1:M]
      residual[1:M] .-= du[1:M]
      residual[1:M] .-= π_tilde # discretized system of ODE for v, where v'(T) = 0 (eq:24)
      residual[M+1] = u[1] + x - dot(ω, u[1:M]) # residual (eq:25)
      residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31)
  end
  @test @inferred f(du, u, p, t) == f(du, u, p, t)
end
