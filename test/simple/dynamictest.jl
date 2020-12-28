
# Setting defaults 
settings = settings_simple(
    z_ex = range(0.0, stop = 5.0, length = 100),
    ts = range(0.0, stop = 100, length = 10) 
) # everything unmentioned is default 

params = parameters_simple(
    υ = 0.1,
    θ = 2.1,
    ζ = 14.5,
    r = 0.05,
    ξ = 1.0
)

params_func = parameters_simple(
    μ = 0.0,
    υ = 0.1, 
    θ = 2.1,
    r = t -> 0.05,
    x = t -> 14.5,
    π = (t, z) -> 1)

params_func_varying_1 = merge(params_func, (π = (t, z) -> 1 + 1e-02 * (1-t)/settings.T,))

params_func_varying_2 = merge(params_func, (r = t -> 0.05 - 1e-02 * (1 - t / settings.T),))

params_func_varying_3 = merge(params_func, (r = t -> 0.05 - 1e-02 * (1 - t / settings.T), π = (t, z) -> 1 + 1e-02 * (1-t)/settings.T))

# Solutions.
# Test the stationary residual is close to zero.
    # even by solving with DAE
    residuals, v_ts, g_ts = solve_simple_transition(params_func, settings)
    @test norm(residuals[1]) ≈ 0 atol = 1e-5
    @test norm(residuals[end]) ≈ 0 atol = 1e-5
    @test norm(residuals) ≈ 0 atol = 1e-5

# Solve with time-varying r and π, now with DAE
    # First experiment. 
    residuals, v_ts, g_ts = solve_simple_transition(params_func_varying_1, settings)
    @test norm(residuals[1]) ≈ 0 atol = 1e-5
    @test norm(residuals[end]) ≈ 0 atol = 1e-5
    @test norm(residuals) ≈ 0 atol = 1e-5

    # Second experiment. 
    residuals, v_ts, g_ts = solve_simple_transition(params_func_varying_2, settings)
    @test norm(residuals[1]) ≈ 0 atol = 1e-5
    @test norm(residuals[end]) ≈ 0 atol = 1e-5
    @test norm(residuals) ≈ 0 atol = 1e-5

    # Third experiment. 
    residuals, v_ts, g_ts = solve_simple_transition(params_func_varying_3, settings)
    @test norm(residuals[1]) ≈ 0 atol = 1e-5
    @test norm(residuals[end]) ≈ 0 atol = 1e-5
    @test norm(residuals) ≈ 0 atol = 1e-5
