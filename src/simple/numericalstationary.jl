# Numerically solve for the stationary solution as a system of equatoins.  
using BandedMatrices, NamedTuples, Roots, Parameters

function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ = params
    
    function stationary_numerical_given_g(g)
        x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
        if γ - g < 0
            L_1 = L_1_minus
        else
            L_1 = L_1_plus
        end
        L_T = (r - g)*I - (γ - g) * L_1 - (σ^2/2.0)* L_2 
        π = exp.(z)
        v_T = L_T \ π
        f_z = α * exp.(-α*z)
        f_i = f_z ./ (1 - exp.(-α*z[end]))
        ω = f_i ./ sum(f_i)


        diff = v_T[1] + ζ - ω' * v_T
        return diff
    end

    g_T = find_zero(stationary_numerical_given_g, (-0.0001, 0.5))

    x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
    if γ - g_T < 0
        L_1 = L_1_minus
    else
        L_1 = L_1_plus
    end
    L_T = (r - g_T)*I - (γ - g_T) * L_1 - (σ^2/2.0) * L_2
    π = exp.(z)
    v_T = L_T \ π
    ν_T = (γ-g_T)/σ^2 + √(((g_T-γ)/σ^2)^2 + (r-g_T)/(σ^2/2));
    return @NT(g = g_T, ν = ν_T, v = v_T)
end 