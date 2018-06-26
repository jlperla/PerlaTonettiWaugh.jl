#Just implements the analytical solution we have in the [notes](https://github.com/jlperla/perla_tonetti_waugh/blob/master/code/dynamics_proof_of_concept/simplified_problem_notes.pdf), returning back functions (g, υ, v) as well as stationary_algebraic. 

# Arnav Sood: Jun. 26, 2018 

using NamedTuples

# Generate params (a function so that it's not top-level).
function simpleparams()
    ans = @NT(γ = 0.005, α = 2.1, ζ = 14.5, r = 0.05, σ = 0.02); # Default values are taken from the `test_stationary_solution.m` file in perla_tonetti_waugh.
    return ans 
end 

# Validate params.
function checkparams(params)
    # Put any parameter restrictions here.
    return true
end

# Validate results
function checkresults(results)
    # Put any results validation here. 
    @assert results.υ >= 0 # Check υ is nonnegative. [EQUATION]
    return true
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
    ans = γ + (1-(α-1)*ζ*(r-γ))/((α-1)^2 * ζ) + σ^2/2 * (α*(α*(α-1)*(r-γ-σ^2/2)*ζ-2)+1)/((α-1)*((α-1)*(r-γ-σ^2/2)*ζ-1)); # Equation 6 
    return ans 
end 

# Calculate stationary υ, given params and g.
function υ(g, params)
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
function v(z, υ, params)
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

# Calculate stationary g and v, given params. 
function stationary_algebraic(params)
    # Validate params. 
    checkparams(params);
    # Calculate g. 
    g_val = g(params); 
    print(g)
    # Calculate υ
    υ_val = υ(g_val, params);
    # Calculate v.
    v_val = z -> v(z, υ_val, params); # Overloads v.
    ans = @NT(g = g_val, υ = υ_val, v = v_val);
    # Check results
    checkresults(ans)
    # Return
    return ans
end 

# Calculate stationary g and v, given params. 
function stationary_algebraic(z, params)
    # Validate params. 
    checkparams(params);
    # Calculate g. 
    g_val = g(params); 
    print(g)
    # Calculate υ
    υ_val = υ(g_val, params);
    # Calculate v.
    v_val = v(z, υ_val, params); # Overloads v.
    ans = @NT(g = g_val, υ = υ_val, v = v_val);
    # Check results
    checkresults(ans)
    # Return 
    return ans
end 
