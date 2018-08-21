function solve_dynamic(params, settings)
    # Setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, T, ode_solve_algorithm, iterations = settings

    g_grid_length = length(z)
    INTERPOLATION_SCHEME_g = LinInterp
    

    # Quadrature weighting
    ω = ω_weights(z, θ, ξ)
    # Time grids
    ts = linspace(0.0, T, g_grid_length)

    # Define and solve dynamic ODE. 
    ode_prob = simpleODE(params, settings)

    # perform transformation to enforce upwind conditions
    g_lb = μ + υ^2/2 + 1e-5
    g_ub = 5 * (μ + υ^2/2)
    transformer = ArrayTransformation(bridge(ℝ, Segment(g_lb, g_ub)), g_grid_length)

    # initial guess for the solution
    g_mid = (g_lb + g_ub) / 2
    g_candidate_initial = inverse(transformer, fill(g_mid, g_grid_length)) 
    
    function dynamic_g_problem_vec!(residuals, g_candidate_vec)
        g_candidate_vec = transformer(g_candidate_vec)
        
        g_candidate_interp = INTERPOLATION_SCHEME_g(ts, g_candidate_vec)
        g_candidate_func = t -> g_candidate_interp(t)
        
        # TODO: there should be a way to do this without unpacking/packing repeatedly like ode_prob.p.g = g_candidate_func
        # at this moment this seems to be the best way to do so as type ODEProblem and one for ode_prob.p are immutable
        @unpack L_1, L_2, z, r, υ, π_tilde, T, μ = ode_prob.p
        p_updated = @NT(L_1 = L_1, L_2 = L_2, z = z, g = g_candidate_func, r = r, υ = υ, π_tilde = π_tilde, T = T, μ = μ)
        ode_prob_updated = ODEProblem(ode_prob.f, ode_prob.u0, (T, 0.0), p_updated)
        
        # compute residuals and update
        residuals_updated = calculate_residuals(ode_prob_updated, x, ω, ode_solve_algorithm, ts)
        residuals[:] = [residual for residual in residuals_updated]
    end
    
    # use nlsolver
    solved = nlsolve(dynamic_g_problem_vec!, g_candidate_initial, iterations = iterations)
    
    residuals = zeros(g_grid_length)
    dynamic_g_problem_vec!(residuals, solved.zero) # compute the resulting ssr
    g = transformer(solved.zero)
    return @NT(g = g, residuals = residuals, solved = solved)
end

