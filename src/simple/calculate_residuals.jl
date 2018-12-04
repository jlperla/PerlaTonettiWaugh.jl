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
    for (i, t) in enumerate(ts)
        v_t = sol(t)[1:M] # i.e., the value function at the point.
        residuals[i] = v_t[1] + x(t) - dot(ω, v_t)
    end
    return residuals
end

function minimize_residuals(params, settings)
    # setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm, t_grid = settings
    g_T = g(T)

    # returns a vector of residuals given a vector of g for linear interpolation
    function calculate_residuals_by_candidate(g_vectorized, params, settings)
        g_interpolated = LinearInterpolation(settings.t_grid, g_vectorized)
        return calculate_residuals(params, merge(settings, (g = g_interpolated, )))
    end

    # setup for optimization
    settings = merge(settings, (t_grid = t_grid, )) # TODO: once passed by settings, remove this line
    g_initial = fill(g_T, length(t_grid))

    # solve the optimization problem
    solved = LeastSquaresOptim.optimize(x -> calculate_residuals_by_candidate(x, params, settings), g_initial, 
                                        autodiff = :central,
                                        LeastSquaresOptim.LevenbergMarquardt(), lower = fill(0.0, length(g_initial)))

    return calculate_residuals_by_candidate(solved.minimizer, params, settings) # return the resultiing residuals 
end