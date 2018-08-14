
# Implementation of the simple model with time-varying objects. 
function simpleODE(params, settings)
    # Unpack necessary objects. 
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, T, g = settings 
    M = length(z) 
    
    # Discretize the operator. 
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, ξ) # L_1_minus ≡ L_1 is the only one we use. 

    # Calculate the stationary solution.
    r_T = r(T)
    g_T = g(T)
    π_tilde_T = z -> π_tilde(T, z)

    # NOTE: the following two lines overlap with stationary solution
    L_T = (r_T - g_T - ξ*((μ+υ^2/2) - g_T) - υ^2/2*ξ^2)*I - ((μ+υ^2/2) - g_T + υ^2*ξ)*L_1_minus - υ^2/2 * L_2 # Construct the aggregate operator.
    v_T = L_T \ π_tilde_T.(z) # Solution to the rescaled differential equation.

    # Bundle as before. 
    p = @NT(L_1 = L_1_minus, L_2 = L_2, z = z, g = g, r = r, υ = υ, π_tilde = π_tilde, T = T, μ = μ) #Named tuple for parameters.

    # Dynamic calculations, defined for each time ∈ t.  
    function f(du,u,p,t)
        @unpack L_1, L_2, z, r, μ, g, υ, π_tilde, T = p 
        # Validate upwind scheme direction. 
        ((μ+υ^2/2) - g(t)) < 0 || error("μ - g must be strictly negative at all times")
        # Carry out calculations. 
        L = (r(t) - g(t) - ξ*((μ+υ^2/2) - g(t)) - υ^2/2*ξ^2)*I - ((μ+υ^2/2) - g(t) + υ^2*ξ)*L_1 - υ^2/2 * L_2
        A_mul_B!(du,L,u)
        du .-= π_tilde.(t, z)
    end

    return ODEProblem(f, v_T, (T, 0.0), p)
end
