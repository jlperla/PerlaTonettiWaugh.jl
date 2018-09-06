# Algebraically solve for the stationary solution.
function stationary_algebraic_simple(params)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ = params
    # Create γ
    γ = μ + υ^2/2
    # Validate parameters.
    # Calculate g.
    g = γ + (1-(θ-1)*ζ*(r-γ))/((θ-1)^2 * ζ) + υ^2/2 * (θ*(θ*(θ-1)*(r-γ-υ^2/2)*ζ-2)+1)/((θ-1)*((θ-1)*(r-γ-υ^2/2)*ζ-1)); # Equation 6
    # Calculate ν
    ν = (γ-g)/υ^2 + √(((g-γ)/υ^2)^2 + (r-g)/(υ^2/2)); # Equation 8
    # Calculate a generic v.
    v(z) = (r - γ - υ^2/2)^(-1) * (exp(z) + 1/ν * exp(-ν*z)); # Equation 7
    # Validate parameters.
    # Return.
    return (g = g, ν = ν, v = v)
end

# Numerically solve for the stationary solution as a system of equations.
function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ, ξ, π_tilde = params
    z, L_1_minus, L_1_plus, L_2  = rescaled_diffusionoperators(z, ξ) #Discretize the operator

    # Define the pdf of the truncated exponential distribution
    ω = ω_weights(z, θ, ξ)

    # Function we're solving. 
    function stationary_numerical_given_g(g)
        L_T = (r - g - ξ*((μ + υ^2/2) - g) - υ^2/2*ξ^2)*I - ((μ + υ^2/2) - g + υ^2*ξ)*L_1_minus - υ^2/2 * L_2 # Construct the aggregate operator.
        v_T = L_T \ π_tilde.(z) # Solution to the rescaled differential equation.
        diff = v_T[1] + ζ - dot(ω, v_T)
        return diff
    end

    # Find the root. 
    g_T = find_zero(stationary_numerical_given_g, (1e-10, 0.75*r), atol = 1e-10, rtol = 1e-10, xatol = 1e-10, xrtol = 1e-10)

    # Check that the solution makes sense.
    @assert((μ + υ^2/2) - g_T < 0)  # Error if γ - g ≡ (μ + υ^2/2) - g_T is positive

    # Recreate what the ODE returned for the value function.
    L_T = (r - g_T - ξ*((μ + υ^2/2) - g_T) - υ^2/2*ξ^2)*I - ((μ + υ^2/2) - g_T + υ^2*ξ)*L_1_minus - υ^2/2 * L_2 # Construct the aggregate operator.
    v_T = L_T \ π_tilde.(z) # Solution to the rescaled differential equation.
    return (g = g_T, v = v_T)
end