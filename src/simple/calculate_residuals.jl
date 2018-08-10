# This function calculates the residual at all points given on g.
function calculate_residuals(params, settings) # To keep the params consistent with other tuples. 
    # Setup
    @unpack γ, σ, α, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm = settings 

    # Quadrature weighting
    ω = ω_weights(z, α, ξ)

    # Define and solve dynamic ODE. 
    ode_prob = simpleODE(params, settings)
    
    return calculate_residuals(ode_prob, x, ω, ode_solve_algorithm)
end

function calculate_residuals(ode_prob, x, ω, ode_solve_algorithm) # To keep the params consistent with other tuples. 
    # Solve ode
    sol = Sundials.solve(ode_prob, ode_solve_algorithm)
    t_vals = sol.t

    # Calculate the residual at each time point
    residuals = zeros(length(t_vals))
    for (i, t) in enumerate(t_vals)
        v_t = sol(t) # i.e., the value function at the point.
        residuals[i] = v_t[1] + x(t) - dot(ω, v_t)
    end
    
    return residuals
end