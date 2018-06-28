# #Implements the algebraic stationary solution from the full model, returning back [WHAT IT RETURNS].

# using Parameters, NamedTuples

# function stationary_algebraic_full(params)
#     # Unpack parameters.
#     @unpack ρ, σ, n, θ, γ, d, κ, ζ, η, Θ, χ, υ, μ, δ = params

#     # Validate parameters.
#     @assert υ > 0 # To avoid discontinuity, by assumption.

#     # Define auxiliary objects.
#     F(z) = 1 - z^(-θ); # (H.1)
#     S = θ * (g - μ - θ * v^2/2); # (H.2)
#     ν = (μ - g)/v^2 + √(((g-μ)/(v^2))^2 + (r-g)/(v^2/2)); # (H.3)
#     a = (r - g - (σ - 1)*(μ - g - (σ-1)v^2/2))^(-1); # (H.4)
#     b = (1 - a*(r-g))*d^(1-σ)*z^(ν + σ - 1); # (H.5)
#     r = ρ + γ*g + δ; # (H.6)
#     L = Ω*((N-1)*(1-F(z))*κ + (1-η)*ζ*(S + δ/χ))
#     # (H.8)
#     ̂z = d * (κ/π)^(1/(σ-1)); # (H.9)
#     w = (̄σ)^(-1) * ̄z; # (H.10)
#     x = ζ * (1 - η + η * Θ/w); # (H.11)
#     # Validate auxiliary objects. 
#     # Define equilibrium equations. 
#     # Solve equilibrium equations. 
#     # Validate equilibrium results.
#     @assert g*(1-γ) < ρ # To ensure finite utility; see Appendix D.  
#     # Return results. 
#     return true
# end 