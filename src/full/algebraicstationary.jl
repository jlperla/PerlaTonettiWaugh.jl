#Implements the algebraic stationary solution from the full model, returning back [WHAT IT RETURNS].

using Parameters, NamedTuples

function stationary_algebraic_full(params)
    # Unpack parameters.
    @unpack ρ, σ, n, θ, γ, d, κ, ζ, η, Θ, χ, υ, μ, δ = params

    # Validate parameters.
    @assert υ > 0 # To avoid discontinuity, by assumption.

    # Define auxiliary objects.
    F(z) = 1 - z^(-θ); # (H.1)
    S(g) = θ * (g - μ - θ * v^2/2); # (H.2)
    r(g) = ρ + γ*g + δ; # (H.6)
    x = ζ * (1 - η + η * Θ/w); # (H.11)
    # Validate auxiliary objects. 
    # Define equilibrium equations. 
    # Solve equilibrium equations. 
    # Validate equilibrium results.
    @assert g*(1-γ) < ρ # To ensure finite utility; see Appendix D.  
    # Return results. 
    return true
end 