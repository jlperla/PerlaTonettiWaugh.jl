# Numerically solve for the stationary solution as a system of equatoins.  
using BandedMatrices, NamedTuples, Roots, Parameters

function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ = params
    function stationary_numerical_given_g(g)
        x, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
        L_T = (r - g)*I - Diagonal(γ - g) * L_1_plus - Diagonal(σ^2/2.0) * L_2 
        π = exp(z)
        v_T = L_T \ π
        f_z = α * exp(-α*z)
        ω = f_z ./ sum(f_z)

        diff = v_T[1] + ζ - ω * v_T
        return diff
    end

f(g) =  stationary_numerical_given_g(g)
g_T = find_zero(f, (0.0001, 2))

@unpack γ, σ, α, r, ζ = params
        x, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
        L_T = (r - g_T)*I - Diagonal(γ - g_T) * L_1_plus - Diagonal(σ^2/2.0) * L_2 
        π = exp(z)
        v_T = L_T \ π
    return v_T
end