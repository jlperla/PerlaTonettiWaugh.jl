# This function calculatesthe residual at all points given on g.
function calculate_residuals(params::NamedTuple, π::Function, g_function::Function, x::AbstractArray, T::Integer) # To keep the params consistent with other tuples. 
    # Setup
    @unpack γ, σ, α, r, ζ = params
    @assert isa(ζ, Function) && isa(r, Function) # Assert that the functional parameters have the right type.
    @assert isa(γ, Number) && isa(σ, Number) && isa(α, Number) # Assert that the constant parameters have the right type. 
    settings=@NT(x = x, g = g_function, T = T, π = π)
    basealgorithm = CVODE_BDF()

    # Quadrature weighting
    ourDist = Truncated(Exponential(1/α), x[1], x[end])
    ω = irregulartrapezoidweights(x, ourDist)

    # Define and solve dynamic ODE. 
    prob = simpledynamicODEproblem(params, settings)
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