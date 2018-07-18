# This function calculatesthe residual at all points given on g.
using NamedTuples, QuantEcon, Parameters

function calculate_residual(params, g_function, z, T)

    # Setup
    @unpack γ, σ, α, r, ζ = params
    settings=@NT(z = z,g = g_function, T = T)
    basealgorithm = CVODE_BDF()

    ζ(t) = ζ

    # Quadrature weighting
    ourDist = Truncated(Exponential(1/α), z[1], z[end])
    ω = irregulartrapezoidweights(z, ourDist)

    prob = create_dynamic_ODE(params, settings)
    sol = solve(prob, basealgorithm)
    t_vals = sol.t

    # calculate the residual at each time point
    for (i, t) in enumerate(t_vals)
        v_t = sol[t]    #i.e., the value function at the point.
        resid[i] = v_t[1] + ζ(t) - dot(ω, v_t)
    end
    
    return resid
end