# Residuals function. 
function entry_residuals(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)
  # Get the solver solution. 
    sol = solve_dynamics(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)
  # Grab the entry residuals and time points. 
    entry_residuals = sol.results[:entry_residual]
    ts = sol.results[:t]
  # Interpolate these over t. 
    entry_residuals_interp = LinInterp(reverse(ts), reverse(entry_residuals)) # This will be returned. 
end 

# Method for interpolation with Ω. 
function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)
  # Validate the interpolation objects.
    @assert Ω_nodes[1] ≈ 0.0 
    @assert Ω_nodes[end] ≈ T
  # Get the (callable) interpolation object. 
    Ω = LinInterp(Ω_nodes, Ω_vec) # QuantEcon routine. 
  # Run the main method. 
    r = solve_dynamics(params_T, stationary_sol_T, settings, T, Ω) 
end 

# Main method.
function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)
    # Unpack arguments 
      @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T # Parameters
      @unpack z, tstops, Δ_E = settings # Settings 
      v_T = stationary_sol_T.v_tilde # Stationary -- 
      g_T = stationary_sol_T.g
      z_hat_T = stationary_sol_T.z_hat
      L_tilde_T = stationary_sol_T.L_tilde
      Ω_T = stationary_sol_T.Ω # -- Stationary 

    # Validate arguments 
      @assert params_T.γ ≈ 1 # γ has to be close 1 to have consistent results with the stationary solutions

    # Define the results data frame we'll be using and push the stationary onto it. 
      results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_0 = v_T[1], L_tilde = L_tilde_T)

    # Define intermediate quantitities. 
      M = length(z)
      ω = ω_weights(z, θ, σ-1) # Quadrature weights. 
      z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operator Discretization. 
      L_1 = L_1_minus # L_1 ≡ L_1_minus.

    # Define E(t) based on forward difference 
      E(t) = (log(Ω(t+Δ_E)) - log(Ω(t-Δ_E)))/(2*Δ_E) + δ

    # Define the intermediate quantities for the DAE problem. 
      S(g) = θ * (g - μ - θ * υ^2/2) # Compute S given g. [EQUATION NUMBER NEEDED]
      L_tilde(S, z_hat, E_t, Ω_t) = Ω_t * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E_t / χ)) # Compute L_tilde. [EQUATION NUMBER NEEDED]

    # Define a function to compute the static equilibrium quantities. 
      function static_equilibrium(v_0, g, z_hat, E_t, Ω_t)
        S_t = S(g)
        L_tilde_t = L_tilde(S_t, z_hat, E_t, Ω_t)
        z_bar = Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)) # [EQUATION NUMBER NEEDED]
        π_min = (1 - L_tilde_t) / ((σ-1)*z_bar) # [EQUATION NUMBER NEEDED]
        i_vectorized = z .>= log(z_hat) # Vectorized indicator function 
        π_tilde = π_min * (1.0.+(N-1)*d^(1-σ)*i_vectorized) - (N-1)*κ*exp.(-(σ-1).*z).*i_vectorized # [EQUATION NUMBER NEEDED]
        entry_residual = v_0 - ζ * (1-χ) / χ # [EQUATION NUMBER NEEDED]
        return (S_t = S_t, L_tilde_t = L_tilde_t, z_bar = z_bar, π_min = π_min, π_tilde = π_tilde, entry_residual = entry_residual)
      end 

    # Set the initial conditions. 
      u0 = [v_T; g_T; z_hat_T]
      du0 = zeros(M+2)

    # Write the actual f! to be called by the DAE solver. 
      function f!(residual,du,u,p,t)
        # Setup (unpack arguments, reset residual, grab E and Ω evaluations, etc.) 
          residual .= 0
          # v = u[1:M] (this line is commented to cut down on memory allocations)
          g = u[M+1]
          z_hat = u[M+2] 
          x = ζ 
          Ω_t = Ω(t)
          E_t = E(t)
        # Get static equilibrium values 
          @unpack S_t, L_tilde_t, z_bar, π_min, π_tilde = static_equilibrium(u[1], g, z_hat, E_t, Ω_t)
        # Grab the L_tilde derivative. 
          L_tilde_log_derivative = 0.0 # Default to Float literal. 
          if (t < T)
            t_forward = results[:t][end]
            L_tilde_forward = results[:L_tilde][end]
            L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde_t))/(t_forward - t) # Forward differences.  
          end 
        #=  Reset the residuals to slack in the DAE conditions. 
            Note that A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2 and we're decomposing this. 
        =#
          residual[1:M] = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:M] # system of ODEs (eq:28)
          residual[1:M] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:M] 
          residual[1:M] .-= (υ^2/2)*L_2*u[1:M] 
          residual[1:M] .-= du[1:M]
          residual[1:M] .-= π_tilde
          residual[M+1] = u[1] + x - dot(ω, u[1:M]) # residual (eq:25)
          residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31)
      end

    # Bundle all of this into an actual DAE problem. 
      dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), differential_vars = [trues(M); false; false])
    
    # Define the callback we'll be using.
      function cb_aux(u, t, integrator) # Function we'll be executing
        # Unpack the u
          # t = t
          g_t = u[M+1]
          z_hat_t = u[M+2]
          Ω_t = Ω(t)
          E_t = E(t)
          v_0_t = u[1]
          # Calculate L_tilde 
            S_t = S(g_t)
            L_tilde_t = L_tilde(S_t, z_hat_t, E_t, Ω_t)
        # Push to results 
          push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_0 = v_0_t, L_tilde = L_tilde_t))
      end 

      cb = FunctionCallingCallback(cb_aux, tdir = -1, func_start = false) # Callback object. 

    # Solve that DAE problem. 
      sol = DifferentialEquations.solve(dae_prob, callback = cb) # solve! # TODO: include tstops = tstops in the argument.

    # Post-process the results DataFrame.
      # Define the welfare, etc. quantities in terms of quantities in the DataFrame. 
        gen_λ_ii = z_hat -> 1 / (1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # [EQUATION NUMBER NEEDED]
        gen_c = (L_tilde, Ω, λ_ii) -> (θ / (1-σ+θ))^(1/(σ-1))*(1-L_tilde)*Ω^(1/(σ-1))*λ_ii^(1/(1-σ)) # [EQUATION NUMBER NEEDED]
        gen_S = S
        gen_z_bar = (Ω_t, z_hat) -> Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)) # [EQUATION NUMBER NEEDED]
        gen_π_min = (L_tilde_t, z_bar) -> (1 - L_tilde_t) / ((σ-1)*z_bar) # [EQUATION NUMBER NEEDED]
        gen_entry_residual = (v_0) -> v_0 - ζ*(1-χ)/χ # [EQUATION NUMBER NEEDED]
      # Add these quantities to the DataFrame. 
        results = @transform(results, λ_ii = gen_λ_ii.(:z_hat)) # λ_ii column. 
        results = @transform(results, c = gen_c.(:L_tilde, :Ω, :z_hat)) # c column.
        results = @transform(results, S = gen_S.(:g)) # S column.
        results = @transform(results, z_bar = gen_z_bar.(:Ω, :z_hat)) # z_bar column. 
        results = @transform(results, π_min = gen_π_min.(:L_tilde, :z_bar)) # π_min column. 
        results = @transform(results, entry_residual = gen_entry_residual.(:v_0)) # entry_residual column 

    # Return. 
      return (results = sort!(results), sol = sol, f! = f!, static_equilibrium = static_equilibrium) # The results, raw DAE solution, and DAE problem (f!, static_equilibrium, etc.) objects.  
end