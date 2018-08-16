function solve_dynamic(params, settings)
    # Setup
    @unpack γ, σ, α, r, x, ξ, π_tilde = params
    @unpack z, T, ode_solve_algorithm = settings

    g_grid_length = length(z)
    INTERPOLATION_SCHEME_g = LinInterp
    INITIAL_CANDIDATES_COUNT = 10 # number of initial guesses

    # Quadrature weighting
    ω = ω_weights(z, α, ξ)
    # Time grids
    ts = linspace(0.0, T, g_grid_length)

    # Define and solve dynamic ODE. 
    ode_prob = simpleODE(params, settings)

    # perform transformation to enforce upwind conditions
    # NOTE: the ub is set to be 5 * γ, but this might not be a good bound 
    g_lb = γ + 1e-5
    g_ub = 5 * γ
    transformer = ArrayTransformation(bridge(ℝ, Segment(g_lb, g_ub)), g_grid_length)

    # generate initial candidates for solutions 
    g_candidate_initials_elements = linspace(g_lb, g_ub, (INITIAL_CANDIDATES_COUNT+2))
    g_candidate_initials_elements = g_candidate_initials_elements[2:(INITIAL_CANDIDATES_COUNT+1)] # discard the most lower & upper candidates
    g_candidate_initials = [fill(v, g_grid_length) for v in g_candidate_initials_elements]
    results_g = zeros(g_grid_length, INITIAL_CANDIDATES_COUNT) # resulting g per each guess
    results_residuals = zeros(g_grid_length, INITIAL_CANDIDATES_COUNT) # resulting residuals per each guess
    ssrs = fill(Inf, INITIAL_CANDIDATES_COUNT) # results are sorted by ssr    
    
    # define optimization problem for g
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
    
    # solve for g from each candidate
    for (i, g_candidate_initial) in enumerate(g_candidate_initials)
        g_candidate_initial = inverse(transformer, g_candidate_initial) 
        solved = solve_dynamic_single(g_candidate_initial, dynamic_g_problem_vec!, transformer)
        ssrs[i] = norm(solved.residuals)
        results_g[:,i] = solved.g
        results_residuals[:,i] = solved.residuals
    end

    # choose the solution with the least ssr
    g = results_g[:,indmin(ssrs)]
    residuals = results_residuals[:,indmin(ssrs)]

    return @NT(g = g, residuals = residuals)
end

function solve_dynamic_single(g_candidate_initial, dynamic_g_problem_vec!, transformer)
    # solve for g from each cnadidate
    solved = nlsolve(dynamic_g_problem_vec!, g_candidate_initial)

    # compute the resulting residuals
    residuals = zeros(length(g_candidate_initial))
    dynamic_g_problem_vec!(residuals, solved.zero) 

    # transform g back
    g = transformer(solved.zero)
    return @NT(g = g, residuals = residuals)
end