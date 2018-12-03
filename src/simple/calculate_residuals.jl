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
