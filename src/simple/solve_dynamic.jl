function solve_dynamic(params, settings)
    # Setup
    @unpack γ, σ, α, r, x, ξ, π_tilde = params
    @unpack z, T, ode_solve_algorithm = settings

    g_grid_length = length(z)
    INTERPOLATION_SCHEME_g = LinInterp
    

    # Quadrature weighting
    ω = ω_weights(z, α, ξ)
    # Time grids
    ts = linspace(0.0, T, g_grid_length)

    # Define and solve dynamic ODE. 
    ode_prob = simpleODE(params, settings)

    # perform transformation to enforce upwind conditions
    # NOTE: the ub is set to be 3*ode_prob.p.g_T, but this might not be a good bound 
    transformer = ArrayTransformation(bridge(ℝ, Segment(γ, (3*ode_prob.p.g_T))), g_grid_length)

    # initial guess for the solution
    g_candidate_initial = inverse(transformer, fill(ode_prob.p.g_T, g_grid_length)) 
    
    function dynamic_g_problem_vec!(residuals, g_candidate_vec)
        g_candidate_vec = transformer(g_candidate_vec)
        
        g_candidate_interp = INTERPOLATION_SCHEME_g(ts, g_candidate_vec)
        g_candidate_func = t -> g_candidate_interp(t)
        
        # TODO: there should be a way to do this without unpacking/packing repeatedly like ode_prob.p.g = g_candidate_func
        # at this moment this seems to be the best way to do so as type ODEProblem and one for ode_prob.p are immutable
        @unpack L_1, L_2, z, r, σ, π_tilde, T, γ, g_T = ode_prob.p
        p_updated = @NT(L_1 = L_1, L_2 = L_2, z = z, g = g_candidate_func, r = r, σ = σ, π_tilde = π_tilde, T = T, γ = γ, g_T = g_T)
        ode_prob_updated = ODEProblem(ode_prob.f, ode_prob.u0, (T, 0.0), p_updated)
        
        # compute residuals and update
        residuals_updated = calculate_residuals(ode_prob_updated, x, ω, ode_solve_algorithm, ts)
        residuals[:] = [residual for residual in residuals_updated]
    end
    
    # use nlsolver
    solved = nlsolve(dynamic_g_problem_vec!, g_candidate_initial)
    
    residuals = zeros(g_grid_length)
    dynamic_g_problem_vec!(residuals, solved.zero) # compute the resulting ssr
    g = transformer(solved.zero)
    return @NT(g = g, residuals = residuals, solved = solved)
end

