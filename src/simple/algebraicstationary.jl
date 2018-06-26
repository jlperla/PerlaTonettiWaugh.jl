#Implements the algebraic stationary solution from the [notes](https://github.com/jlperla/perla_tonetti_waugh/blob/master/code/dynamics_proof_of_concept/simplified_problem_notes.pdf), returning back functions (g, υ, v) as well as stationary_algebraic. 

using Parameters, NamedTuples

function stationary_algebraic_simple(params)
    # Unpack parameters. 
    @unpack γ, σ, α, r, ζ = params
    # Validate parameters.
    # Calculate g. 
    g = γ + (1-(α-1)*ζ*(r-γ))/((α-1)^2 * ζ) + σ^2/2 * (α*(α*(α-1)*(r-γ-σ^2/2)*ζ-2)+1)/((α-1)*((α-1)*(r-γ-σ^2/2)*ζ-1)); # Equation 6 
    # Calculate υ
    υ = (γ-g)/σ^2 + √(((g-γ)/σ^2)^2 + (r-g)/(σ^2/2)); # Equation 8
    # Calculate a generic v. 
    v(z) = (r - γ - σ^2/2)^(-1) * (exp(z) + 1/υ * exp(-υ*z)); # Equation 7
    # Validate parameters. 
    # Return. 
    return @NT(g = g, υ = υ, v = v)
end