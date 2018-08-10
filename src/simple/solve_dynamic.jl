function solve_dynamic(params, settings)
    # Setup
    @unpack γ, σ, α, r, x, ξ, π_tilde = params
    @unpack z, T, ode_solve_algorithm = settings
    # TODO: let users choose the number of grids in settings param later
    GRID_LENGTH_g = 1000
    # TODO: let users choose interpolation scheme
    INTERPOLATION_SCHEME = LinInterp

    # Quadrature weighting
    ω = ω_weights(z, α, ξ)
    # Time grids
    ts = linspace(0.0, T, GRID_LENGTH_g)

    # Define and solve dynamic ODE. 
    ode_prob = simpleODE(params, settings)

    function dynamic_g_problem(g_candidate_vec)
        g_candidate_interp = INTERPOLATION_SCHEME(ts, g_candidate_vec)
        g_candidate_func = t -> g_candidate_interp(t)

        # TODO: there should be a way to do this without unpacking/packing repeatedly like ode_prob.p.g = g_candidate_func
        # at this moment this seems to be the best way to do so as type ODEProblem and one for ode_prob.p are immutable
        @unpack L_1, L_2, z, r, σ, π_tilde, T, γ, g_T = ode_prob.p
        p_updated = @NT(L_1 = L_1, L_2 = L_2, z = z, g = g_candidate_func, r = r, σ = σ, π_tilde = π_tilde, T = T, γ = γ, g_T = g_T)
        ode_prob_updated = ODEProblem(ode_prob.f, ode_prob.u0, (T, 0.0), p_updated)
        
        residuals = calculate_residuals(ode_prob_updated, x, ω, ode_solve_algorithm)
        
        ssr = 0
        # this is faster than norm(residuals)
        for i = 1:length(residuals)
            ssr += residuals[i] * residuals[i]
        end

        return ssr 
    end
    
    g_candidate_initial = convert(Array{Float64}, linspace(ode_prob.p.g_T, ode_prob.p.g_T, GRID_LENGTH_g))
    
    # FIXIT: this does not work due to constraints for upwind scheme -- use JuMP.jl?
    solved = optimize(dynamic_g_problem, g_candidate_initial, BFGS())
    
    g_equilibrium_vec = solved.minimizer
    g_equilibrium_interp = INTERPOLATION_SCHEME(ts, g_equilibrium_vec)
    g_equilibrium = t -> g_equilibrium_interp(t)

    return @NT(g = g_equilibrium, ssr = solved.minimum)
end