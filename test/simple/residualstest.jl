# User settings.
# State grid.
    z_min = 0.0
    z_max = 5.0
    M = 100
    z_grid = range(z_min, stop = z_max, length = M) # Since we only care about the grid.

# Time grid.
    T_val = 100.0
    N = 10
    t = range(0.0, stop = T_val, length = N) # For interpolation purposes only.

# Constant parameters.
    υ_val = 0.1
    θ_val = 2.1
    ζ_val = 14.5
    r_val = 0.05
    γ_val = 0.005
    ξ_val = 1.0

# Functional parameters.
    x_func = t -> ζ_val # Not idiosyncratic, per equation (4)
    r_func_varying = t -> (r_val - 1e-02 * (1 - t / T_val)) # Not idiosyncratic, per intro to doc.
    r_func_ = t -> r_val

# π_tilde_func_varying = (t, z) -> r_val .- 1e-02 * (1 .- 0* t / T_val) # Not idiosyncratic, per intro to doc.
    π_tilde_func_varying = (t, z) -> (1 + 1e-02 * (1 - t / T_val))
    π_tilde_func_ = (t, z) -> 1 # Potentially idiosyncratic.

# Param generators and param NTs.
    params_ = @with_kw (μ = 0.0, υ = υ_val, θ = θ_val, r = r_val, ζ = ζ_val, ξ = ξ_val, π_tilde = z -> 1)
    params_func = @with_kw (μ = 0.0, υ = υ_val, θ = θ_val, r = r_func_, x = x_func, ξ = ξ_val, π_tilde = π_tilde_func_)
    params_ = params_()
    params_func_ = params_func()
    params_func_varying_1 = params_func(π_tilde = π_tilde_func_varying)
    params_func_varying_2 = params_func(r = r_func_varying)
    params_func_varying_3 = params_func(r = r_func_varying, π_tilde = π_tilde_func_varying)

# Solutions.
# Solve for the numerical stationary g_T.
    result_ns = stationary_numerical_simple(params_, z_grid)
    g_stationary = result_ns.g # This is the level.

# Create the interpolation object of g
    g_vector = g_stationary .+ 0.01 * t
    g_int = LinearInterpolation(t, g_vector)
    g_func = t -> g_int(t) # Not idiosyncratic.

# Create settings object.
    settings = @with_kw (z = z_grid, T = T_val, g = t -> g_stationary, 
                        ode_solve_algorithm = CVODE_BDF(), iterations = 1000, 
                        t_grid = range(0.0, stop = T_val, length = length(z_grid)),
                        g_node_count = 30)

# Test the stationary residual is close to zero.
    resid = calculate_residuals(params_func_, settings())
    @test norm(resid) ≈ 0 atol = 1e-5
    # even by solving with DAE
    ω = ω_weights(z_grid, θ_val, ξ_val)
    daeprob = simpleDAE(params_func_, settings())
    resid = calculate_residuals(daeprob, x_func, ω, IDA(), t)
    @test norm(resid[1]) ≈ 0 atol = 1e-5
    @test norm(resid[end]) ≈ 0 atol = 1e-5
    @test norm(resid) ≈ 0 atol = 1e-5

# Solve with time-varying r and π_tilde, now with DAE
    daeprob = simpleDAE(params_func_varying_1, settings())
    resid = calculate_residuals(daeprob, x_func, ω, IDA(), t)
    @test norm(resid[1]) ≈ 0 atol = 1e-5
    @test norm(resid[end]) ≈ 0 atol = 1e-5
    @test norm(resid) ≈ 0 atol = 1e-5
    daeprob = simpleDAE(params_func_varying_2, settings())
    resid = calculate_residuals(daeprob, x_func, ω, IDA(), t)
    @test norm(resid[1]) ≈ 0 atol = 1e-5
    @test norm(resid[end]) ≈ 0 atol = 1e-5
    @test norm(resid) ≈ 0 atol = 1e-5
    daeprob = simpleDAE(params_func_varying_3, settings())
    resid = calculate_residuals(daeprob, x_func, ω, IDA(), t)
    @test norm(resid[1]) ≈ 0 atol = 1e-5
    @test norm(resid[end]) ≈ 0 atol = 1e-5
    @test norm(resid) ≈ 0 atol = 1e-5

@testset "minimize_residuals by least SSR" begin 
    # try minimize_residuals by least square optimization instead of solving ODE by DAE
    # need higher tolerance (atol = 1e-3)
    resid = minimize_residuals(params_func_varying_1, settings())
    @test norm(resid[1]) ≈ 0 atol = 1e-3
    @test norm(resid[end]) ≈ 0 atol = 1e-3
    @test norm(resid) ≈ 0 atol = 1e-3
    resid = minimize_residuals(params_func_varying_2, settings())
    @test norm(resid[1]) ≈ 0 atol = 1e-3
    @test norm(resid[end]) ≈ 0 atol = 1e-3
    @test norm(resid) ≈ 0 atol = 1e-3
    resid = minimize_residuals(params_func_varying_3, settings())
    @test norm(resid[1]) ≈ 0 atol = 1e-3
    @test norm(resid[end]) ≈ 0 atol = 1e-3
    @test norm(resid) ≈ 0 atol = 1e-3
end