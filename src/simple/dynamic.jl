
# Implementation of the simple model with time-varying objects. 
function simpleODE(params, settings)
    # Unpack necessary objects. 
    @unpack γ, σ, ζ, r, α = params
    @unpack π, x, T, g = settings 
    M = length(x)

    # Discretize the operator. 
    x, L_1, L_1_plus, L_2 = diffusionoperators(x) # L_1_minus ≡ L_1 is the only one we use. 

    #Calculate the stationary solution.
    (γ - g(T)) < 0 || error("γ - g must be strictly negative at all times") # How we validate the upwind scheme. 
    L_T = (r(T) - g(T))*I - (γ - g(T)) * L_1 - σ^2/2.0 * L_2
    v_T = L_T \ π.(T, x)
    @assert(issorted(v_T)) # Monotonicity check. 

    # Bundle as before. 
    p = @NT(L_1 = L_1, L_2 = L_2, x = x, g = g, r = r, σ = σ, π = π, T = T, γ = γ) #Named tuple for parameters.

    # Dynamic calculations, defined for each time ∈ t.  
    function f(du,u,p,t)
        @unpack L_1, L_2, x, r, γ, g, σ, π, T = p 
        # Validate upwind scheme direction. 
        (γ - g(t)) < 0 || error("μ - g must be strictly negative at all times")
        # Carry out calculations. 
        L = (r(t) - g(t))*I - (γ - g(t)) * L_1 - σ^2/2.0 * L_2 # Aggregated discrete operator. 
        A_mul_B!(du,L,u)
        du .-= π.(t, x)
    end

    #Checks on the stationary residual
    dv_T = zeros(v_T)
    f(dv_T, v_T, p, T)
    @show norm(dv_T)
    @assert norm(dv_T) < 1.0E-10

    return ODEProblem(f, v_T, (T, 0.0), p)
end
