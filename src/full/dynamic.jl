# Kernel function for main method.
  function f!(residual,du,u,p,t)
    # Setup (unpack arguments, reset residual, grab E and Ω evaluations, etc.)
      @unpack ζ, Ω, E, static_equilibrium, T, results, ρ, δ, σ, μ, υ, L_1, L_2, ω, κ, d = p
      residual .= 0
      P = length(residual) - 2
      g = u[P+1]
      z_hat = u[P+2]
      x = ζ
      Ω_t = Ω(t)
      E_t = E(t)
    # Get static equilibrium values
      @unpack S_t, L_tilde_t, z_bar, π_min, π = static_equilibrium(u[1], g, z_hat, E_t, Ω_t)
    # Grab the L_tilde derivative.
      L_tilde_log_derivative = 0.0 # Default to Float literal.
      if (t < T)
        t_forward = results[:t][end]
        L_tilde_forward = results[:L_tilde][end]
        L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde_t))/(t_forward - t) # See note under (34)
      end
    #=  Reset the residuals to slack in the DAE conditions.
        Note that (C.40) and (46) yield A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2 and we're decomposing this.
    =#
      residual[1:P] = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:P] # (34)
      residual[1:P] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:P] # (46)
      residual[1:P] .-= (υ^2/2)*L_2*u[1:P] # (46)
      residual[1:P] .-= du[1:P]
      residual[1:P] .-= π # discretized system of ODE for v, where v'(T) = 0 (47)
      residual[P+1] = u[1] + x - dot(ω, u[1:P]) # value matching residual, (48) and x(t) = ζ assumption at beginning of Section 2
      residual[P+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (49)
  end

# Main method.
function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; detailed_solution = true)

    if (T < settings.T_U_bar)
      throw("Terminal time `T` should be large enough so that T >= settings.T_U_bar is satisfied.")
    end

    # Unpack arguments
      @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T # Parameters
      @unpack z, T_U_bar, tstops = settings # Settings
      v_T = stationary_sol_T.v_tilde # Stationary --
      g_T = stationary_sol_T.g
      z_hat_T = stationary_sol_T.z_hat
      L_tilde_T = stationary_sol_T.L_tilde
      Ω_T = stationary_sol_T.Ω # -- Stationary

    # Validate arguments
      @assert γ ≈ 1 # γ has to be close 1 to have consistent results with the stationary solutions
      @assert η == 0

    # Define the results data frame we'll be using and push the stationary onto it.
      results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_0 = v_T[1], L_tilde = L_tilde_T)

    # Define intermediate quantitities.
      P = length(z)
      ω = ω_weights(z, θ, σ-1) # Quadrature weights.
      z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operator Discretization.
      L_1 = L_1_minus # L_1 ≡ L_1_minus.

    # Define the auxiliary functions for the DAE problem.
      S(g) = θ * (g - μ - θ * υ^2/2) # Compute S given g. (26)
      L_tilde(S, z_hat, E_t, Ω_t) = Ω_t * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E_t / χ)) # Compute L_tilde. (27)

      function static_equilibrium(v_0, g, z_hat, E_t, Ω_t)
        S_t = S(g)
        L_tilde_t = L_tilde(S_t, z_hat, E_t, Ω_t)
        z_bar = Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)) # (31)
        π_min = (1 - L_tilde_t) / ((σ-1)*z_bar) # (32)
        i_vectorized = z .>= log(z_hat) # Vectorized indicator function
        π = π_min * (1.0.+(N-1)*d^(1-σ)*i_vectorized) - (N-1)*κ*exp.(-(σ-1).*z).*i_vectorized # (33)
        entry_residual = v_0 - ζ * (1-χ) / χ # value matching condition (45)
        return (S_t = S_t, L_tilde_t = L_tilde_t, z_bar = z_bar, π_min = π_min, π = π, entry_residual = entry_residual)
      end

    # Set the initial conditions.
      u0 = [v_T; g_T; z_hat_T]
      du0 = zeros(P+2)

    # Create the parameters object
      p = (ζ = ζ, Ω = Ω, E = E, static_equilibrium = static_equilibrium, T = T,
            results = results, ρ = ρ, δ = δ, σ = σ, μ = μ, υ = υ, L_1 = L_1, L_2 = L_2,
            ω = ω, κ = κ, d = d)

    # Bundle all of this into an actual DAE problem.
      dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false])

    # Define the callback we'll be using (i.e., for the backward-looking L_tilde derivative)
      function cb_aux(u, t, integrator) # Function we'll be executing
        # Unpack the u
          # t = t
          g_t = u[P+1]
          z_hat_t = u[P+2]
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
      if (tstops == nothing)
        sol = DifferentialEquations.solve(dae_prob, callback = cb)
      else
        sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops)
      end

    # Post-process the results DataFrame.
    results = sort!(results)
      # Define the welfare, etc. quantities in terms of quantities in the DataFrame.
        gen_λ_ii = z_hat -> 1 / (1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (51)
        gen_c = (L_tilde, Ω, z_bar, S) -> (1 - L_tilde)*z_bar - η*ζ*Ω*Theta*(S + δ/χ) # (52)
        gen_S = S
        gen_z_bar = (Ω_t, z_hat) -> (Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)))^(1/(σ-1)) # (31)
        gen_π_min = (L_tilde_t, z_bar) -> (1 - L_tilde_t) / ((σ-1)*z_bar) # (32)
        gen_entry_residual = (v_0) -> v_0 - ζ*(1-χ)/χ # (45)
        gen_L_tilde_adopt = (Ω, S) -> Ω * ζ * S # (30)
        gen_L_tilde_export = (Ω, z_hat) -> Ω * ((N-1)*z_hat^(-θ))*κ # (28)
        gen_L_tilde_entrycost = (Ω, E) -> Ω * ζ * E / χ # (29)

      # Add these quantities to the DataFrame.
        results = @transform(results, entry_residual = gen_entry_residual.(:v_0)) # entry_residual column

        log_c_T = log(gen_c(L_tilde_T, Ω_T, gen_z_bar(Ω_T, z_hat_T), S(g_T)))

        g_interpolated(t) = (sol(t))[P+1]
        z_hat_interpolated(t) = (sol(t))[end]
        L_tilde_interpolated(t) = L_tilde(S(g_interpolated(t)), z_hat_interpolated(t), E(t), Ω(t))
        z_bar(t) = gen_z_bar(Ω(t), z_hat_interpolated(t))
        U(t) = quadgk(τ -> exp(-ρ*τ)*(log_M(t+τ) + log_c(t+τ)), 0, (T-t))[1] + exp(-ρ*(T-t))*(g_T + ρ*(log_c_T + g_T * T))/(ρ^2)
        c(t) = gen_c(L_tilde_interpolated(t), Ω(t), z_bar(t), S(t))
        log_M(t) = quadgk(g_interpolated, 0, t)[1]
        log_c(t) = log(gen_c(L_tilde_interpolated(t), Ω(t), gen_z_bar(Ω(t), z_hat_interpolated(t)), S(g_interpolated(t))))

        U_bar_T_generator(t, T_cutoff) = quadgk(τ -> exp(-ρ*τ)*(log_M(t+τ) + log_c(t+τ)), 0, (T_cutoff-t))[1] + exp(-ρ*(T_cutoff-t))*(g_T + ρ*(log_c_T + g_T * T_cutoff))/(ρ^2)
        U_bar_T(t) = U_bar_T_generator(t, T_U_bar)
        if (detailed_solution)
          # other welfare functions.
          # L_tilde_interpolated = LinearInterpolation(results[:t], results[:L_tilde])
          λ_ii = t -> gen_λ_ii(z_hat_interpolated(t))

          results = @transform(results, λ_ii = gen_λ_ii.(:z_hat)) # λ_ii column.
          results = @transform(results, S = gen_S.(:g)) # S column.
          results = @transform(results, z_bar = gen_z_bar.(:Ω, :z_hat)) # z_bar column.
          results = @transform(results, c = gen_c.(:L_tilde, :Ω, :z_bar, :S)) # c column.
          results = @transform(results, π_min = gen_π_min.(:L_tilde, :z_bar)) # π_min column.
          results = @transform(results, log_M = log_M.(:t)) # log_M column
          results = @transform(results, U = U.(:t)) # U column
          results = @transform(results, L_tilde_a = gen_L_tilde_adopt.(:Ω, :S))
          results = @transform(results, L_tilde_x = gen_L_tilde_export.(:Ω, :z_hat))
          results = @transform(results, L_tilde_E = gen_L_tilde_entrycost.(:Ω, :E))
        end

    # Return.
    # The results, raw DAE solution, and DAE problem (f!, static_equilibrium, etc.) objects.
      return (results = results, sol = sol, p = p, static_equilibrium = static_equilibrium, 
              U = U, c = c, Ω = Ω, log_M = log_M, log_c = log_c, U_bar_T = U_bar_T,
              U_bar_T_generator = U_bar_T_generator) 
end
 