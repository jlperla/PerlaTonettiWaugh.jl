# Global constants. 
    # State grid. 
    x_min = 0.0 # Rename to X for consistency. 
    x_max = 5.0
    M = 100
    x_grid = linspace(x_min,x_max,M) # Since we only care about the grid. 

    # Time grid. 
    T_val = 10.0
    N = 10
    t = linspace(0.0, T_val, N) # For interpolation purposes only. 

    # Functional parameters. 
    π_func = (t, x) -> exp.(x) # Potentially idiosyncratic. 
    ζ_func = t -> ζ_val # Not idiosyncratic, per equation (4)
    r_func = t -> r_val # Not idiosyncratic, per intro to doc. 

    # Constant parameters. 
    σ_val = 0.02
    α_val = 2.1
    ζ_val = 14.5
    r_val = 0.05
    γ_val = 0.005
    ξ_val = 1.0
    
    # Param generators and param NTs. 
    params = @with_kw (γ = γ_val, r = r_val, ζ = ζ_val, α = α_val, σ = σ_val, ξ = ξ_val, π_tilde = x -> 1) # Callable generator. 
    params_const = params()
    params_func = params(r = r_func, ζ = ζ_func) 

# Solutions.
    # Solve for the numerical stationary g_T. 
    result_ns = stationary_numerical_simple(params_const, x_grid)
    g_stationary = result_ns.g # This is the level. 

    # Create the interpolation object of g
    g_vector = g_stationary + 0.01 * t
    g_int = LinInterp(t, g_vector)
    g_func = t -> g_int(t) # Not idiosyncratic. 

    # Create settings object.
    settings = @with_kw (x = x_grid, T = T_val, π = π_func, g = g_func)

    # Calculate residuals. 
    resid = calculate_residuals(params_func, settings())
    @show resid # Need this one for now, I think. 

    # Test the stationary residual is close to zero.
    resid2 = calculate_residuals(params_func, settings(g = t -> g_stationary))
    @test_broken norm(resid2) ≈ 0 atol = 1e-10