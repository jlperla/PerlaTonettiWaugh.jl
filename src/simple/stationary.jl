# Algebraically solve for the stationary solution.
function stationary_algebraic_simple(params)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ = params
    # Calculate g.
    g = μ + (1-(θ-1)*ζ*(r-μ))/((θ-1)^2 * ζ) + υ^2/2 * (θ*(θ*(θ-1)*(r-μ-υ^2/2)*ζ-2)+1)/((θ-1)*((θ-1)*(r-μ-υ^2/2)*ζ-1)); # (9)
    # Calculate ν
    ν = (μ-g)/υ^2 + √(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)); # (11)
    # Calculate a generic v.
    v(z) = 1/(r - μ - υ^2/2) * (1 + 1/ν * exp(-(ν+1)*z)); # (10). Note that these values are in the transformed space.
    # Validate parameters.
    # Return.
    return (g = g, ν = ν, v = v)
end

# Numerically solve for the stationary solution as a system of equations.
function stationary_numerical_simple(params, z_ex)
    @assert z_ex[1] == 0.0 && issorted(z_ex)  # validate grid
    z = z_ex[2:end-1]

    # Unpack parameters
    @unpack μ, υ, θ, r, ζ, π = params

    # Define the pdf of the truncated exponential distribution
    ω = ω_weights(z_ex, θ, 1) # (26), see utils/quadrature.jl

    # Differential objects
    bc = (Mixed(1), Mixed(1)) # boundary conditions for differential operators
    L_1_minus = L₁₋(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂(z_ex, bc)

    Ξ₁ = 1/(1 - 1*(z[1] - 0.0)) # (24)

    # Function we're solving
    function stationary_numerical_given_g(in)
        g = in[1]
        # Construct the aggregate operator.
        A = (r - μ - υ^2/2)*I - (μ + υ^2 - g)*L_1_minus - υ^2/2 * L_2 # (19)
        v = A \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (30)
        diff = Ξ₁*v[1] - dot(ω, v) + ζ # value matching condition (31)
        return diff
    end

    # Find and validate the root
    sol = solve_system(stationary_numerical_given_g, [0.1])
    g_T = sol[1]
    @assert(μ + υ^2/2 - g_T < 0) # Negative drift condition (20)

    # Use the g_T to recreate L_T and v_T
    A_T = (r - μ - υ^2/2)*I - (μ + υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (19)
    v_T = A_T \ π.(z) # (30)
    return (g = g_T, v = v_T)
end
