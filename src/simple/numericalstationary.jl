# Numerically solve for the stationary solution as a system of equatoins.  
using BandedMatrices, NamedTuples, Roots, Parameters

function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ = params

    # Define the pdf of the truncated exponential distribution
    ourDist = Truncated(Exponential(α), z[1], z[end])
    ω = irregulartrapezoidweights(z, ourDist)

    function stationary_numerical_given_g(g)
        x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator

        L_T = (r - g)*I - (γ - g) * L_1_minus - (σ^2/2.0)* L_2 
        π = exp.(z)
        v_T = L_T \ π

        diff = v_T[1] + ζ - ω' * v_T
        return diff
    end

    g_T = find_zero(stationary_numerical_given_g, (1e-10, 0.75*r))

    #assert(γ - g_T < 0)  # Error if γ - g is positive

    x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
    L_T = (r - g_T)*I - (γ - g_T) * L_1_minus - (σ^2/2.0) * L_2
    π = exp.(z)
    v_T = L_T \ π
    ν_T = (γ-g_T)/σ^2 + √(((g_T-γ)/σ^2)^2 + (r-g_T)/(σ^2/2));
    return @NT(g = g_T, ν = ν_T, v = v_T)
end 