# This function calculates the residual at all points given on g.
function calculate_residuals(params, settings) # To keep the params consistent with other tuples.
    # Setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm, t_grid = settings
    # Calculations
    ω = ω_weights(z, θ, ξ) # Quadrature weighting
    ode_prob = simpleODE(params, settings)  # Define and solve dynamic ODE.
    return calculate_residuals(ode_prob, x, ω, ode_solve_algorithm, t_grid)
end

function calculate_residuals(ode_prob, x, ω, ode_solve_algorithm, ts) # To keep the params consistent with other tuples.
    # Solve ode
    sol = Sundials.solve(ode_prob, ode_solve_algorithm, tstops = ts)
    M = length(ω)
    # Calculate the residual at each time point
    residuals = zeros(length(ts))
    v_ts = zeros(M, length(ts))
    for (i, t) in enumerate(ts)
        v_t = sol(t)[1:M] # i.e., the value function at the point.
        residuals[i] = v_t[1] + x(t) - dot(ω, v_t) # (eq:A.20)
        v_ts[:,i] = v_t
    end

    return (residuals = residuals, v_ts = v_ts)
end

function minimize_residuals(params, settings)
    # setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm, t_grid, g_node_count = settings
    g_T = g(T)
    interpolate_g(g_vectorized, settings) = CubicSplineInterpolation(settings.t_node_for_g, [g_vectorized; settings.g_T])
    # returns a vector of residuals given a vector of g for linear interpolation
    function calculate_residuals_by_candidate(g_vectorized, params, settings)
        residuals, v_ts = calculate_residuals(params, merge(settings, (g = interpolate_g(g_vectorized, settings), )))
        return residuals
    end
    # setup for optimization
    settings = merge(settings, (t_node_for_g = range(0.0, stop = T, length = g_node_count), g_T = g_T))
    g_initial = fill(g_T, g_node_count - 1)
    # solve the optimization problem
    solved = LeastSquaresOptim.optimize(x -> calculate_residuals_by_candidate(x, params, settings), g_initial, autodiff = :central, LeastSquaresOptim.LevenbergMarquardt(), lower = fill(0.0, length(g_initial)))
    settings = merge(settings, (g = interpolate_g(solved.minimizer, settings), ))
    return calculate_residuals(params, settings) # return the resultiing residuals
end

function minimize_residuals_python(params, settings)
    # setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm, t_grid, g_node_count = settings
    g_T = g(T)
    interpolate_g(g_vectorized, settings) = CubicSplineInterpolation(settings.t_node_for_g, [g_vectorized; settings.g_T])
    # returns a vector of residuals given a vector of g for linear interpolation
    function calculate_residuals_by_candidate(g_vectorized, params, settings)
        residuals, v_ts = calculate_residuals(params, merge(settings, (g = interpolate_g(g_vectorized, settings), )))
        return residuals
    end
    # setup for optimization
    settings = merge(settings, (t_node_for_g = range(0.0, stop = T, length = g_node_count), g_T = g_T))
    g_initial = fill(g_T, g_node_count - 1)
    # solve the optimization problem
    sol = DFOLS.solve(x -> calculate_residuals_by_candidate(x, params, settings), g_initial, bounds = (fill(0., size(g_initial)), nothing))
    settings = merge(settings, (g = interpolate_g(DFOLS.optimizer(sol), settings), )) # only runs if converged(sol), otherwise throws error(sol.msg)
    return calculate_residuals(params, settings) # return the resultiing residuals
end
