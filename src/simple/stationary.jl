function stationary_algebraic_simple(parameters, settings)
    @unpack μ, υ, θ, r, ζ = parameters
    @unpack T = settings
    g = μ + (1 - (θ - 1)*ζ*(r(T) - μ))/((θ - 1)^2 * ζ) + υ^2/2 * (θ*(θ*(θ - 1)*(r(T) - μ - υ^2/2)*ζ - 2) + 1)/((θ - 1)*((θ - 1)*(r(T) - μ - υ^2/2)*ζ - 1)); # (9)
    ν = (μ-g)/υ^2 + √(((g-μ)/υ^2)^2 + (r(T) - g)/(υ^2/2)); # (11)
    v(z) = 1/(r(T) - μ - υ^2/2) * (1 + 1/ν * exp(-(ν + 1)*z)); # (10). Note that these values are in the transformed space.
    return (g = g, ν = ν, v = v)
end

# Minimizes value-matching condition (31)
function stationary_numerical_simple(parameters, settings)
    @unpack z_ex, T = settings
    x0 = settings.stationary_x0(parameters, settings)
    @unpack μ, υ, θ, r, ζ, π = parameters
    @assert z_ex[1] == 0.0 && issorted(z_ex)  # validate grid
    z = z_ex[2:end-1]
    r_T = r(T)

    ω = ω_weights(z_ex, θ, 1) # (26), see utils/quadrature.jl

    bc = (Mixed(1), Mixed(1)) # boundary conditions for differential operators
    L_1_minus = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)

    Ξ₁ = 1/(1 - 1*(z[1] - 0.0)) # (24)

    function stationary_numerical_given_g(in)
        g = in[1]
        A = (r_T - μ - υ^2/2)*I - (μ + υ^2 - g)*L_1_minus - υ^2/2 * L_2 # (19)
        v = A \ π.(T, z) # discretized system of ODE for v, where v'(T) = 0 (30)
        return Ξ₁*v[1] - dot(ω, v) + ζ # value matching condition (31)
    end

    sol = nlsolve(stationary_numerical_given_g, x0, inplace = false)
    g_T = sol.zero[1]
    @assert converged(sol) "The solver didn't converge; exiting now."
    @assert (μ + υ^2/2 - g_T < 0) # Negative drift condition (20)
    @assert r_T > g_T "r_T > g_T is required for convergence" # From (19)

    # Use the g_T to recreate L_T and v_T
    A_T = (r_T - μ - υ^2/2)*I - (μ + υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (19)
    v_T = A_T \ π.(T, z) # (30)
    return (g = g_T, v = v_T, residual = stationary_numerical_given_g(g_T))
end
