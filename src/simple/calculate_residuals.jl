# This function calculatesthe residual at all points given on g.
function calculate_residuals(params, settings) # To keep the params consistent with other tuples. 
    # Setup
    @unpack γ, σ, α, r, ζ, ξ, π_tilde = params
    @unpack x, g, T = settings 
    # Solver setting 
    basealgorithm = CVODE_BDF()

    # Quadrature weighting
    ω = ω_weights(x, α, ξ)

    # Define and solve dynamic ODE. 
    prob = simpleODE(params, settings)
    sol = Sundials.solve(prob, basealgorithm)
    t_vals = sol.t

    # Calculate the residual at each time point
    resid = zeros(length(t_vals))
    for (i, t) in enumerate(t_vals)
        v_t = sol(t)    #i.e., the value function at the point.
        resid[i] = v_t[1] + ζ(t) - dot(ω, v_t)
    end
    
    return resid
end