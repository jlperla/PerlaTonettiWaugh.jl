# This function calculatesthe residual at all points given on g.
using NamedTuples, QuantEcon, Parameters

function calculate_residuals(params, g_function, z, T)
    # Setup
    @unpack γ, σ, α, r, ζ, π = params # We don't need pi explicitly here, but the create_dynamic_ODE() call requires it. 
    settings=@NT(z = z, g = g_function, T = T)
    basealgorithm = CVODE_BDF()

    # Quadrature weighting
    ourDist = Truncated(Exponential(1/α), z[1], z[end])
    ω = irregulartrapezoidweights(z, ourDist)

    prob = create_dynamic_ODE(params, settings)
    sol = Sundials.solve(prob, basealgorithm)
    t_vals = sol.t

    # calculate the residual at each time point
    resid = zeros(length(t_vals))
    for (i, t) in enumerate(t_vals)
        v_t = sol(t)    #i.e., the value function at the point.
        resid[i] = v_t[1] + ζ(t) - dot(ω, v_t)
    end
    
    return resid
end