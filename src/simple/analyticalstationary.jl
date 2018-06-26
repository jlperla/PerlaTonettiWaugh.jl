#Just implements the analytical solution we have in the notes returning back functions and g

using NamedTuples
params = @NT(γ = 2, α = 2, ζ = 2, r = 2, σ = 2)

# Validate params.
function checkparams(params)
    # Put any parameter restrictions here.
end 

# Calculate stationary g, given params. 
function g(params)
    # Unpack params. 
    γ = params.γ; 
    α = params.α; 
    ζ = params.ζ; 
    r = params.r; 
    σ = params.σ;
    # Calculate g. 
    ans = γ + (1-(α-1)*ζ*(r-γ))/((α-1)^2 * ζ) + σ^2/2 * (α(α*(α-1)*(r-γ-σ^2/2)*ζ-2)+1)/((α-1)*((α-1)*(r-γ-σ^2/2)*ζ-1)); # Equation 6 
    return ans 
end 

# Calculate stationary υ, given params and g.
function υ(params, g)
     # Unpack params. 
     γ = params.γ; 
     α = params.α; 
     ζ = params.ζ; 
     r = params.r; 
     σ = params.σ;
     # Calculate υ
     ans = (γ-g)/σ^2 + √(((g-γ)/σ^2)^2 + (r-g)/(σ^2/2)); # Equation 8
     return ans
end

# Calculate v, given z, params and υ. 
function v(z, params, g, υ)
     # Unpack params. 
     γ = params.γ; 
     α = params.α; 
     ζ = params.ζ; 
     r = params.r; 
     σ = params.σ;
     # Calculate v 
     ans = (r - γ - σ^2/2)^(-1) * (exp(z) + 1/υ * exp(-υ*z)); # Equation 7
     return ans 
end

# Calculate stationary g and v, given params and a state. 
function stationary_algebraic(params, z)
    # Validate params. 
    checkparams(params);
    # Calculate g. 
    g = g(params); # Overloads g. 
    # Calculate υ
    υ = υ(params, g); # Overloads υ. 
    # Calculate v.
    v = v(z, params, g, υ); # Overloads v.
    ans = @NT(g = g, υ = υ, v = v);
    return ans
end 

