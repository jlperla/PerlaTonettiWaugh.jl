# This function calculatesthe residual at all points given on g.
function calculate_residuals(params, settings) # To keep the params consistent with other tuples. 
    # Setup
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, g, T, ode_solve_algorithm = settings 

    # Quadrature weighting
    ω = ω_weights(z, θ, ξ)

    # Define and solve dynamic ODE. 
    prob = simpleODE(params, settings)
    sol = Sundials.solve(prob, ode_solve_algorithm)
    t_vals = sol.t

    # Calculate the residual at each time point
    residuals = zeros(length(t_vals))
    for (i, t) in enumerate(t_vals)
        v_t = sol(t)    #i.e., the value function at the point.
        residuals[i] = v_t[1] + x(t) - dot(ω, v_t)
    end
    
    return residuals
end