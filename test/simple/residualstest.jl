# Global constants. 
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
    settings = @with_kw (z = z_grid, T = T_val, g = t -> g_stationary, ode_solve_algorithm = CVODE_BDF(), iterations = 1000)
    # Solve for v with time-varying g
    resid = calculate_residuals(params_func_, settings(g = g_func))

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
  
#=
    Test advanced options.
=#

    # tstops. 
    daeprob = simpleDAE(params_func_, settings()) # Use the vanilla example.
    sol_vanilla = DifferentialEquations.solve(daeprob)
    tstops = [1.0, 11.3, 27.4, 30.0] # Arbitrary tstops. 
    sol_tstops = DifferentialEquations.solve(daeprob, tstops = tstops)
    @test tstops ⊆ sol_tstops.t # Check to see if the tstops are included. 
    @test norm(sol_vanilla.u[end] - sol_tstops.u[end]) ≈ 0.0 atol = 1e-5 # Check that the two solutions end up fairly close to one another. 

    # differencing.
        # Naive way to do it (no callbacks).
        function f(u, p, t) # Default example from DifferentialEquations.jl tutorial. 
            val = 1.01*u
            push!(p, (t, val))
            return val
        end 
        u0=1/2
        tspan = (0.0,-1.0)
        tstops = [-0.1, -0.2, -0.3, -1.0]
        p = []
        odeprob = ODEProblem(f,u0,tspan,p)
        sol_vanilla = DifferentialEquations.solve(odeprob, tstops = tstops)

        function growth(p, tstops)
            ourPoints = [p[i] for i = 1:2:length(p) if p[i][1] ∈ tstops] # Need the 2 because the solutions seem to come in slightly different pairs. 
            diffs = [(ourPoints[i+1][2] - ourPoints[i][2])/(ourPoints[i+1][1] - ourPoints[i][1]) for i = 1:length(ourPoints)-1]
        end 

        # proper way to do it (with callbacks).
        u0=1/2
        tspan = (0.0,-1.0)
        tstops = [-0.1, -0.2, -0.3, -1.0]
        p = (saved_values = SavedValues(Float64, Tuple{Float64,Float64}),)
        odeprob = ODEProblem((u, p, t) -> begin
            vals = p.saved_values.saveval
                if t < 0.0
                i = findlast(x -> x[1] > t, vals); # Sign is becuase of backwards time. 
                d = (u - vals[i][2])/(vals[i][1] - t); # Differences for backwards time, per notes. 
            end 
            return 1.01*u 
        end, u0, tspan, p) # To confirm that we have access to saved_values during the runs. 
        cb = SavingCallback((u,t,integrator)->(t, u), p.saved_values, tdir = -1)
        sol_callback = DifferentialEquations.solve(odeprob, callback=cb, tstops = tstops)
        print(p.saved_values.saveval)
        print(p.saved_values.t)



