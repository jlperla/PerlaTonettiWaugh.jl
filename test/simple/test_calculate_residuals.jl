# Global constants. 
    # State grid. 
    x_min = 0.0 # Rename to X for consistency. 
    x_max = 5.0
    M = 100
    x=linspace(x_min,x_max,M) # Since we only care about the grid. 
    # Time grid. 
    T = 10.0
    N = 10
    t = linspace(0.0, T, N)
    # Functions, params, etc. 
    π_func(t, x) = exp(x)*(1.0+0.0*t);
    σ_n = 0.02
    α_val = 2.1
    ζ_n = 14.5
    r_n = 0.05
    ζ_func(t) = ζ_n + 0.0*t
    r_func(t, x) = r_n + 0.0*x + 0.0*t
    σ_new(t, x) = σ_n + 0.0*x + 0.0*t
    γ=0.005
    # Param wrappers. 
    params_dynamic = @with_kw (γ=0.005, σ=σ_new, π=π_func, ζ=ζ_func, r = r_func, α=α_val)
    params_constant = @with_kw (γ=0.005, σ=0.02, α=2.1, r=0.05, ζ=14.5)

# Solutions. 
    # solve for numerical g_T as test
    res = stationary_numerical_simple(params_constant(),z)
    g_analytic=res.g # This, I think, is just wrong. The result returns the numerical g? 
    v_analytic=res.v
    @show v_analytic

    # Create the interpolation object of g
    g_vector = g_analytic + 0.01 * t
    g_int=LinInterp(t, g_vector)
    g(t, x) = g_int(t) + 0.0 * x

    # Calculate residuals. 
    resid = calculate_residuals(params_dynamic(), g, z, T)
    @show resid