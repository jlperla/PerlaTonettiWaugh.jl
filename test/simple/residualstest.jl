# Global constants. 
    # State grid. 
    z_min = 0.0 
    z_max = 5.0
    M = 100
    z_grid = linspace(z_min,z_max,M) # Since we only care about the grid. 

    # Time grid. 
    T_val = 100.0
    N = 10
    t = linspace(0.0, T_val, N) # For interpolation purposes only. 

    # Functional parameters. 
    π_tilde_func = (t, z) -> 1 # Potentially idiosyncratic. 
    x_func = t -> ζ_val # Not idiosyncratic, per equation (4)
    r_func = t -> r_val # Not idiosyncratic, per intro to doc. 

    # Constant parameters. 
    σ_val = 0.02
    α_val = 2.1
    ζ_val = 14.5
    r_val = 0.05
    γ_val = 0.005
    ξ_val = 1.0
    
    # Param generators and param NTs. 
    params_const = @with_kw (γ = γ_val, σ = σ_val, α = α_val, r = r_val, ζ = ζ_val, ξ = ξ_val, π_tilde = z -> 1) 
    params_func = @with_kw (γ = γ_val, σ = σ_val, α = α_val, r = r_func, x = x_func, ξ = ξ_val, π_tilde = π_tilde_func)
    params_const = params_const()
    params_func = params_func()
 
    # Solutions.
    # Solve for the numerical stationary g_T. 
    result_ns = stationary_numerical_simple(params_const, z_grid)
    g_stationary = result_ns.g # This is the level. 

    # Create the interpolation object of g
    g_vector = g_stationary + 0.01 * t
    g_int = LinInterp(t, g_vector)
    g_func = t -> g_int(t) # Not idiosyncratic. 

    # Create settings object.
    settings = @with_kw (z = z_grid, T = T_val, g = g_func, ode_solve_algorithm = CVODE_BDF())

    # Solve for v with time-varying g
    resid1 = calculate_residuals(params_func, settings())
    @test_broken norm(resid) ≈ 0 atol = 1e-10 # since time-varying g is not in equilibrium, we expect broken at this moment

    # Test the stationary residual is close to zero.
    resid2 = calculate_residuals(params_func, settings(g = t -> g_stationary))
    @test norm(resid2) ≈ 0 atol = 1e-10