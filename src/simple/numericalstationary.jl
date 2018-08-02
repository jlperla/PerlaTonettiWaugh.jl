# Numerically solve for the stationary solution as a system of equatoins.  
function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ, ξ = params

    # Define the pdf of the truncated exponential distribution
    ourDist = Truncated(Exponential(1/α), z[1], z[end])
    ω = irregulartrapezoidweights(z, ourDist)

    function stationary_numerical_given_g(g)
        z, L_1_minus, L_1_plus, L_2  = rescaled_diffusionoperators(z, ξ) #Discretize the operator

        r_tild = r - g - ξ*(γ - g) - σ^2/2*ξ^2;
        μ_tild = γ - g + σ^2*ξ;
        π_tild(z) = 1;
        L_T = r_tild*I - μ_tild*L_1_minus - σ^2/2 * L_2 # Construct the aggregate operator. 
        v_T = L_T \ π_tild(z) # Solution to the rescaled differential equation. 

        diff = v_T[1] + ζ - dot(ω, exp.(ξ.*z).*v_T)
        return diff
    end

    g_T = find_zero(stationary_numerical_given_g, (1e-10, 0.75*r), atol = 1e-10, rtol = 1e-10, xatol = 1e-10, xrtol = 1e-10)

    # Check that the solution makes sense. 
    @assert(γ - g_T < 0)  # Error if γ - g is positive

    # Recreate what the ODE returned for the value function. 
    z, L_1_minus, L_1_plus, L_2  = rescaled_diffusionoperators(z, ξ) #Discretize the operator
    r_tild = r - g_T - ξ*(γ - g_T) - σ^2/2*ξ^2;
    μ_tild = γ - g_T + σ^2*ξ;
    π_tild(z) = 1;
    L_T = r_tild*I - μ_tild*L_1_minus - σ^2/2 * L_2 # Construct the aggregate operator. 
    v_T = L_T \ π_tild(z) # Solution to the rescaled differential equation. 
    return @NT(g = g_T, v = v_T)
end 