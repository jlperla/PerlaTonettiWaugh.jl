# Algebraically solve for the stationary solution.
function stationary_algebraic_simple(params)
    # Unpack parameters.
    @unpack γ, σ, α, r, ζ = params
    # Validate parameters.
    # Calculate g.
    g = γ + (1-(α-1)*ζ*(r-γ))/((α-1)^2 * ζ) + σ^2/2 * (α*(α*(α-1)*(r-γ-σ^2/2)*ζ-2)+1)/((α-1)*((α-1)*(r-γ-σ^2/2)*ζ-1)); # Equation 6
    # Calculate ν
    ν = (γ-g)/σ^2 + √(((g-γ)/σ^2)^2 + (r-g)/(σ^2/2)); # Equation 8
    # Calculate a generic v.
    v(z) = (r - γ - σ^2/2)^(-1) * (exp(z) + 1/ν * exp(-ν*z)); # Equation 7
    # Validate parameters.
    # Return.
    return @NT(g = g, ν = ν, v = v)
end

# Numerically solve for the stationary solution as a system of equations.
function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters.
    @unpack γ, σ, α, r, ζ, ξ, π_tilde = params
    x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, ξ) #Discretize the operator

    # Define the pdf of the truncated exponential distribution
    ω = ω_weights(z, α, ξ)

    function stationary_numerical_given_g(g)
        L_T = (r - g - ξ*(γ - g) - σ^2/2*ξ^2)*I - (γ - g + σ^2*ξ)*L_1_minus - σ^2/2 * L_2 # Construct the aggregate operator.
        v_T = L_T \ π_tilde.(z) # Solution to the rescaled differential equation.
        diff = v_T[1] + ζ - dot(ω, v_T)
        return diff
    end

    g_T = find_zero(stationary_numerical_given_g, (1e-10, 0.75*r), atol = 1e-10, rtol = 1e-10, xatol = 1e-10, xrtol = 1e-10)

    # Check that the solution makes sense.
    @assert(γ - g_T < 0)  # Error if γ - g is positive

    # Recreate what the ODE returned for the value function.
    z, L_1_minus, L_1_plus, L_2  = rescaled_diffusionoperators(z, ξ) #Discretize the operator
    L_T = (r - g_T - ξ*(γ - g_T) - σ^2/2*ξ^2)*I - (γ - g_T + σ^2*ξ)*L_1_minus - σ^2/2 * L_2 # Construct the aggregate operator.
    v_T = L_T \ π_tilde.(z) # Solution to the rescaled differential equation.
    return @NT(g = g_T, v = v_T)
end