
# Implementation of the simple model with time-varying objects. 
function simpleODE(params, settings)
    # Unpack necessary objects. 
    @unpack γ, σ, α, r, ζ, ξ, π_tilde = params
    @unpack x, T, g = settings 
    M = length(x)

    # Discretize the operator. 
    x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(x, ξ) # L_1_minus ≡ L_1 is the only one we use. 

    # Calculate the stationary solution.
    r_T = r(T)
    g_T = g(T)
    π_tilde_T = x -> π_tilde(T, x)

    # NOTE: the following two lines overlap with stationary solution
    L_T = (r_T - g_T - ξ*(γ - g_T) - σ^2/2*ξ^2)*I - (γ - g_T + σ^2*ξ)*L_1_minus - σ^2/2 * L_2 # Construct the aggregate operator.
    v_T = L_T \ π_tilde_T.(x) # Solution to the rescaled differential equation.

    # Bundle as before. 
    p = @NT(L_1 = L_1_minus, L_2 = L_2, x = x, g = g, r = r, σ = σ, π_tilde = π_tilde, T = T, γ = γ) #Named tuple for parameters.

    # Dynamic calculations, defined for each time ∈ t.  
    function f(du,u,p,t)
        @unpack L_1, L_2, x, r, γ, g, σ, π_tilde, T = p 
        # Validate upwind scheme direction. 
        (γ - g(t)) < 0 || error("μ - g must be strictly negative at all times")
        # Carry out calculations. 
        L = (r(t) - g(t) - ξ*(γ - g(t)) - σ^2/2*ξ^2)*I - (γ - g(t) + σ^2*ξ)*L_1 - σ^2/2 * L_2
        A_mul_B!(du,L,u)
        du .-= π_tilde.(t, x)
    end

    #Checks on the stationary residual
    dv_T = zeros(v_T)
    f(dv_T, v_T, p, T)
    # @show norm(dv_T)
    # @assert norm(dv_T) < 1.0E-10

    return ODEProblem(f, v_T, (T, 0.0), p)
end
