# Numerically solve for the stationary solution as a system of equatoins.  
using BandedMatrices, NamedTuples, Roots, Parameters

function stationary_numerical_simple(params, z)
    M = length(z)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ = params


    d = diff(z)
    Δ_m = zeros(M)
    Δ_m[1] = d[1]    # using the first difference as diff from ghost node
    Δ_m[2:end] = d
    Δ_p = zeros(M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    ω_bar = (Δ_p+Δ_m)/2
    ω_bar[1] = Δ_p[1]/2
    ω_bar[end] = Δ_m[end]/2
    f_z = α * exp.(-α*z)
    f_i = f_z ./ (1 - exp.(-α*z[end]))
    ω = ω_bar .* f_i

    
    function stationary_numerical_given_g(g)
        x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator

        L_T = (r - g)*I - (γ - g) * L_1_minus - (σ^2/2.0)* L_2 
        π = exp.(z)
        v_T = L_T \ π

        diff = v_T[1] + ζ - ω' * v_T
        return diff
    end

    g_T = find_zero(stationary_numerical_given_g, (-1.0E-10, 0.5))

    #assert(γ - g_T < 0)  # Error if γ - g is positive

    x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(z, M) #Discretize the operator
    L_T = (r - g_T)*I - (γ - g_T) * L_1_minus - (σ^2/2.0) * L_2
    π = exp.(z)
    v_T = L_T \ π
    ν_T = (γ-g_T)/σ^2 + √(((g_T-γ)/σ^2)^2 + (r-g_T)/(σ^2/2));
    return @NT(g = g_T, ν = ν_T, v = v_T)
end 