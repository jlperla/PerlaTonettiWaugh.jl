#= 
    DAE Algorithms.
=#

    # DAE with time-varying objects. 
    function simpleDAEproblem(params::NamedTuple, x::AbstractArray)
        # Unpack parameters. 
        @unpack ρ, c̃, μ̃, σ, T = params

        # Discretize operator.
        x, L_1, L_1_plus, L_2 = diffusionoperators(x) #Discretize the operator
        M = length(x)

        # Solve the stationary problem. 
        μ̃(T) < 0 || error("μ̃ must be strictly negative at all times.") # Check on the upwind direction. 
        L_T = ρ*I - Diagonal(mu_tilde.(T, x)) * L_1 - Diagonal(σ̃.(T, x).^2/2.0) * L_2
        u_T = L_T \ c̃.(T, x)
        u_ex_T = [u_T; 1.0] #Closed form for the trivial linear function we are adding
        @assert(issorted(u_T)) #We are only solving versions that are increasing for now

        # Bundle. 
        p = @NT(L_1 = L_1, L_2 = L_2, x = x, ρ = ρ, μ̃ = μ̃, σ̃ = σ̃, c̃ = c̃, M = M) #Named tuple for parameters.


        function f(resid,du_ex,u_ex,p,t)
            @unpack L_1, L_2, x, ρ, μ̃, σ̃, c̃, M = p 
            all(μ̃.(t, x) .< 0) || error("μ̃ must be strictly negative at all times.") # Check on the upwind direction. 
            L = ρ*I - Diagonal(μ̃.(t, x)) * L_1 - Diagonal(σ̃.(t, x).^2/2.0) * L_2 
            u = u_ex[1:M]
            resid[1:M] .= L * u_ex[1:p.M] - c̃.(t, p.x)
            resid[M+1] .= u_ex[M+1] - 1.0
            resid .-= du_ex
        end

        # Test the stationary residual 
        resid_T = zeros(u_ex_T) 
        du_ex_T = zeros(u_ex_T)
        f(resid_T, du_ex_T, u_ex_T, p, T)
        @show norm(resid_T)
        @assert norm(resid_T) < 1.0E-10

        tspan = (T, 0.0)
        return DAEProblem(f, zeros(u_ex_T), u_ex_T, tspan, differential_vars = [trues(u_T); false], p)
    end

#= 
    ODE Algorithms. 
=# 

    # ODE with constant parameters.
    function simpleODEproblem(params::NamedTuple, x::AbstractArray)
        # Unpack params. 
        @unpack ρ, c̃, μ̃, σ̃, T = params 

        # Discretize operator. 
        x, L_1, L_1_plus, L_2 = diffusionoperators(x)
        
        #Calculate the stationary solution.
        all(μ̃.(T, x) .< 0) || error("μ - g must be strictly negative at all times and states") # How we validate the upwind scheme. 
        L_T = ρ(T)*I - μ̃(T) * L_1 - σ̃^2/2.0 * L_2
        v_T = L_T \ π.(T, x)
        @assert(issorted(v_T)) # Monotonicity check. 
        
        # Bundle as before. 
        p = @NT(L_1 = L_1, L_2 = L_2, x = x, ρ = ρ, μ̃ = μ̃, σ = σ, c̃ = c̃, T = T) #Named tuple for parameters.

        function f(du,u,p,t)
            @unpack L_1, L_2, x, ρ, μ̃, σ, c̃, T = p # Unpack params. 
            all(μ̃.(T, x) .< 0) || error("μ̃ must be strictly negative at all times and states")
            L = ρ*I - Diagonal(μ̃.(t, x)) * L_1 - Diagonal(sigma_tilde.(t, x).^2/2.0) * L_2
            A_mul_B!(du,L,u)
            du .-= c̃.(t, x)
        end

        #Checks on the stationary residual. 
        du_T = zeros(u_T)
        f(du_T, u_T, p, T)
        #@show norm(du_T)
        @assert norm(du_T) < 1.0E-10

        return ODEProblem(f, u_T, (T, 0.0), p)
    end

    # ODE with time-varying objects. 
    function simpledynamicODEproblem(params, settings)
        # Unpack necessary objects. 
        @unpack γ, σ, ζ, r, α = params
        @unpack π, x, T, g = settings 
        M = length(x)

        # Discretize the operator. 
        x, L_1, L_1_plus, L_2 = diffusionoperators(x) # L_1_minus ≡ L_1 is the only one we use. 
        
        # Define mu and rho functions 
        μ̃ = t -> γ - g(t)
        ρ = t -> r(t) - g(t)

        #Calculate the stationary solution.
        μ̃(T) < 0 || error("μ - g must be strictly negative at all times") # How we validate the upwind scheme. 
        L_T = ρ(T)*I - μ̃(T) * L_1 - σ^2/2.0 * L_2
        v_T = L_T \ π.(T, x)
        @assert(issorted(v_T)) # Monotonicity check. 

        # Bundle as before. 
        p = @NT(L_1 = L_1, L_2 = L_2, x = x, ρ = ρ, μ̃ = μ̃, σ = σ, π = π, T = T) #Named tuple for parameters.
    
        # Dynamic calculations, defined for each time ∈ t.  
        function f(du,u,p,t)
            @unpack L_1, L_2, x, ρ, μ̃, σ, π, T = p 
            # Validate upwind scheme direction. 
            μ̃(t) < 0 || error("μ - g must be strictly negative at all times")
            # Carry out calculations. 
            L = ρ(t)*I - μ̃(t) * L_1 - σ^2/2.0 * L_2 # Aggregated discrete operator. 
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
