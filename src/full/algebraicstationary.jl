#Implements the algebraic stationary solution from the full model, returning back [WHAT IT RETURNS].

using NamedTuples

# Generate params. 
function fullparams()
    ans = @NT(ρ = 0.05, σ = 3, n = 10, θ = 3.22, γ = 0, d = 5.49, κ = 0.06, ζ = 1.9, η = 0, Θ = 1, χ = 1/3, υ = 0.01, μ = 0, δ = 0.01); # Default values taken from the `../code/transition_dynamics/default_transition_parameters.m` file in the MATLAB repo. 
    return ans 
end 

# Validate params. 
function checkparams(params)
    # Put any parameter restrictions here. 
    return true
end 

# Validate results. 
function checkresults(results)
    # Put any results validation here. 
    return true 
end

# Calculate g, given params. 
function g(params)
    return true 
end 

# Calculate auxiliary objects. 
function findauxiliaries(params)
    # Unpack params. 
    ρ = params.ρ; 
    σ = params.σ; 
    n = params.n;
    θ = params.θ;  
    γ = params.γ; 
    d = params.d; 
    κ = params.κ; 
    ζ = params.ζ; 
    η = params.η; 
    Θ = params.Θ;
    χ = params.χ; 
    υ = params.υ; 
    μ = params.μ; 
    δ = params.δ; 
    # Calculate objects. 
    F(z) = 1 - z^(-θ) # (H.1)
    S = θ*(g - μ)
end 