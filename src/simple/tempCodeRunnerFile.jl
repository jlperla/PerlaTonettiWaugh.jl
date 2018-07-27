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