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
function stationary_numerical_simple(params, z)
    P = length(z)
    # Unpack parameters.
    @unpack μ, υ, θ, r, ζ, π = params

    bc = (Mixed(1), Mixed(1)) # boundary conditions for differential operators
    z_extended = [z[1] - diff(z)[1]; z; z[end] + diff(z)[end]]
    L_1_minus = L₁₋(z_extended, bc) # use backward difference as the drift is negative
    L_2 = L₂(z_extended, bc) 

    # Define the pdf of the truncated exponential distribution
    ω = ω_weights(z, θ, 1) # (20)
    # Function we're solving.
    function stationary_numerical_given_g(in)
        g = in[1]
        # Construct the aggregate operator.
        A = (r - μ - υ^2/2)*I - (μ + υ^2 - g)*L_1_minus - υ^2/2 * L_2 # (17)
        v = A \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (24)
        diff = v[1] + ζ - dot(ω, v) # value matching condition (23)
        return diff
    end
    # Find and validate the root.
    sol = solve_system(stationary_numerical_given_g, [0.1])
    g_T = sol[1]
    @assert(μ + υ^2/2 - g_T < 0) # Negative drift condition (18)
    # Use the g_T to recreate L_T and v_T.
    A_T = (r - μ - υ^2/2)*I - (μ + υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (17)
    v_T = A_T \ π.(z) # (24)
    return (g = g_T, v = v_T)
end
