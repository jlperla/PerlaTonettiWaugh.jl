# Algebraically solve for the stationary solution.
function stationary_algebraic_simple(params)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ = params
    # Calculate g.
    g = μ + (1-(θ-1)*ζ*(r-μ))/((θ-1)^2 * ζ) + υ^2/2 * (θ*(θ*(θ-1)*(r-μ-υ^2/2)*ζ-2)+1)/((θ-1)*((θ-1)*(r-μ-υ^2/2)*ζ-1)); # (eq:8)
    # Calculate ν
    ν = (μ-g)/υ^2 + √(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)); # (eq:10)
    # Calculate a generic v.
    v(z) = 1/(r - μ - υ^2/2) * (1 + 1/ν * exp(-(ν+1)*z)); # (eq:9). Note that these values are in the transformed space.
    # Validate parameters.
    # Return.
    return (g = g, ν = ν, v = v)
end

# Numerically solve for the stationary solution as a system of equations.
function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ, ξ, π_tilde = params
    z, L_1_minus, L_1_plus, L_2  = rescaled_diffusionoperators(z, ξ) # Discretize the operator
    # Define the pdf of the truncated exponential distribution
    ω = ω_weights(z, θ, ξ)
    # Function we're solving.
    function stationary_numerical_given_g(in)
        g = in[1]
        # Construct the aggregate operator.
        L = (r - g - ξ*(μ - g) - ξ^2 * υ^2/2)*I - (μ + ξ*υ^2 - g)*L_1_minus - (υ^2/2)*L_2 # (eq:A.9)
        v = L \ π_tilde.(z) # discretized system of ODE for v, where v'(T) = 0
        diff = v[1] + ζ - dot(ω, v) # value matching condition (eq:A.20)
        return diff
    end
    # Find and validate the root.
    sol = find_zero(stationary_numerical_given_g, [0.1])
    g_T = sol[1]
    @assert(μ + υ^2/2 - g_T < 0) # Negative drift condition.
    # Use the g_T to recreate L_T and v_T.
    L_T = (r - g_T - ξ*(μ - g_T) - ξ^2 * υ^2/2)*I - (μ + ξ*υ^2 - g_T)*L_1_minus - υ^2/2 * L_2
    v_T = L_T \ π_tilde.(z)
    return (g = g_T, v = v_T)
end
