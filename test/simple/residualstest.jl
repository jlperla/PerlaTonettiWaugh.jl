using PerlaTonettiWaugh, Base.Test
using Distributions, Sundials, BenchmarkTools, QuantEcon, Interpolations, Parameters, NamedTuples, NLsolve, ContinuousTransformations

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

    # Constant parameters. 
    υ_val = 0.1
    θ_val = 2.1
    ζ_val = 14.5
    r_val = 0.05
    γ_val = 0.005
    ξ_val = 1.0 
    
    # Functional parameters. 
    x_func = t -> ζ_val # Not idiosyncratic, per equation (4)

    r_vector = r_val + 1e-04 * (1 - t / T_val)
    r_int = LinInterp(t, r_vector)
    r_func_varying = t -> r_int(t) # Not idiosyncratic, per intro to doc.    
    r_func_const = t -> r_val

    π_tilde_vector = 1 + 1e-04 * (1 - t / T_val)
    π_tilde_int = LinInterp(t, π_tilde_vector)
    π_tilde_func_varying = (t, z) -> π_tilde_int(t) # Not idiosyncratic, per intro to doc.    
    π_tilde_func_const = (t, z) -> 1 # Potentially idiosyncratic. 
    
    # Param generators and param NTs. 
    params_const = @with_kw (μ = 0.0, υ = υ_val, θ = θ_val, r = r_val, ζ = ζ_val, ξ = ξ_val, π_tilde = z -> 1) 
    params_func = @with_kw (μ = 0.0, υ = υ_val, θ = θ_val, r = r_func_const, x = x_func, ξ = ξ_val, π_tilde = π_tilde_func_const)
    params_const = params_const()
    params_func_const = params_func()
    params_func_varying_1 = params_func(π_tilde = π_tilde_func_varying)
    params_func_varying_2 = params_func(r = r_func_varying)
    params_func_varying_3 = params_func(r = r_func_varying, π_tilde = π_tilde_func_varying)

    
    # Solutions.
    # Solve for the numerical stationary g_T. 
    result_ns = stationary_numerical_simple(params_const, z_grid)
    g_stationary = result_ns.g # This is the level. 

    # Create the interpolation object of g
    g_vector = g_stationary + 0.01 * t
    g_int = LinInterp(t, g_vector)
    g_func = t -> g_int(t) # Not idiosyncratic. 

    # Create settings object.
    settings = @with_kw (z = z_grid, T = T_val, g = t -> g_stationary, ode_solve_algorithm = CVODE_BDF())
    # Solve for v with time-varying g
    resid = calculate_residuals(params_func_const, settings(g = g_func))
    @test_broken norm(resid) ≈ 0 atol = 1e-10 # since time-varying g is not in equilibrium, we expect broken at this moment

    # Test the stationary residual is close to zero.
    resid = calculate_residuals(params_func_const, settings())
    @test norm(resid) ≈ 0 atol = 1e-5

    # Test the stationary residual is close to zero when using solve_dynamic
    solved = solve_dynamic(params_func_const, settings())
    @test norm(solved.residuals) ≈ 0 atol = 1e+2

    # Solve with time-varying r and π_tilde
    solved = solve_dynamic(params_func_varying_1, settings())
    @test norm(solved.residuals) ≈ 0 atol = 1e+2
    solved = solve_dynamic(params_func_varying_2, settings())
    @test norm(solved.residuals) ≈ 0 atol = 1e+2
    solved = solve_dynamic(params_func_varying_3, settings())
    @test norm(solved.residuals) ≈ 0 atol = 1e+2