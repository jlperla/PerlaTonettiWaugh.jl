function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)
    # Unpack arguments 
      @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T # Parameters
      @unpack z, tstops, Δ_E = settings # Settings 
      v_T = stationary_sol_T.v_tilde # Stationary -- 
      g_T = stationary_sol_T.g
      z_hat_T = stationary_sol_T.z_hat
      L_tilde_T = stationary_sol_T.L_tilde
      Ω_T = stationary_sol_T.Ω # -- Stationary 

    # Validate parameters 
      @assert params_T.γ ≈ 1 # γ has to be close 1 to have consistent results with the stationary solutions

    # Define intermediate quantitities. 
      M = length(z)
      ω = ω_weights(z, θ, σ-1) # [TODO: compute integral for eq 30] <-- may be an obsolete comment
      z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operator Discretization. 

    # Create the p bundle. 
        p = (L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ,
            ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, ω = ω, Ω = Ω,
            v_T = v_T, g_T = g_T, z_hat_T = z_hat_T, Ω_T = Ω_T)

    # Define E(t) based on forward difference 
      E(t) = (log(Ω(t+Δ_E)) - log(Ω(t-Δ_E)))/(2*Δ_E) + δ

    # Define the corresponding DAE problem 
      dae = PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)

    # Define the results data frame we'll be using and push the stationary onto it. 
      resultsDF = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_0 = v_T[1], L_tilde = L_tilde_T, entry_residual = 0.0)

    # Define the callback we'll be using.
      function cb_aux(u, p, t, integrator) # Function we'll be executing
        # Unpack the u
          # t = t
          g_t = u[M+1]
          z_hat_t = u[M+2]
          Ω_t = p.Ω(t)
          E_t = E(t)
          v_0_t = u[1]
          # Calculate L_tilde 
            S_t = dae.S(g)
            L_tilde_t = dae.L_tilde(S_t, z_hat, E, Ω)
        # Push to results 
          push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, L_tilde = L_tilde_t))
      end 

      cb = FunctionCallingCallback(cb_aux, tdir = -1, func_start = false) # Callback object. 

    # Solve that DAE problem and grab the solutions. 
      sol = DifferentialEquations.solve(dae.dae_prob, callback = cb) # solve! # TODO: include tstops = tstops in the argument.
      
    # Post-process the results dataframe.

    # Return.

end

# Implementation of the full model with time-varying objects, represented by DAE
function PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)
    # Unpack params and settings.
    @unpack z = settings
    M = length(z)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, θ, δ, χ, N, ζ, ρ = p

    get_S(g) = θ * (g - μ - θ * υ^2/2)
    get_L_tilde(S, z_hat, E, Ω) = Ω * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E / χ))

    function stationary_equilibrium(v_1, g, z_hat, E, Ω, t)
        S = get_S(g)
        L_tilde = get_L_tilde(S, z_hat, E, Ω)
        z_bar = Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        π_tilde = π_min * (1.0.+(N-1)*d^(1-σ)*(z .>= log(z_hat))) - (N-1)*κ*exp.(-(σ-1).*z).*(z .>= log(z_hat))
        entry_residual = v_1 - ζ * (1-χ) / χ
        return (S = S, L_tilde = L_tilde, z_bar = z_bar,
                π_min = π_min, π_tilde = π_tilde,
                entry_residual = entry_residual)
    end

    # Dynamic calculations, defined for each time ∈ t.
    function f!(residual,du,u,p,t)
        residual[:] = zeros(M+2)

        # Carry out calculations.
        # v = u[1:M]; note that this line is disabled and `v` is replaced with u[1:M] to prevent huge memory allocation 
        g = u[M+1]
        z_hat = u[M+2]

        x = ζ
        Ω_t = Ω(t)
        E_t = E(t)

        S = get_S(g)
        L_tilde = get_L_tilde(S, z_hat, E_t, Ω_t)
        z_bar = Ω(t) * (θ / (1.0 + θ - σ)) * (1.0 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        π_tilde = π_min * (1.0.+(N-1)*d^(1-σ)*(z .>= log(z_hat))) - (N-1)*κ*exp.(-(σ-1).*z).*(z .>= log(z_hat))

        # compute the derivative of L_tilde
        L_tilde_derivative_term = 0
        if (t < T)
            t_forward = p.ts[end]
            L_tilde_forward = p.vals[end]
            L_tilde_derivative_term = (log(1 - L_tilde_forward) - log(1-L_tilde))/(t_forward-t) # Reverse direction. 
        end

        # Form the DAE at t.
        # Note that
        # A_t = (ρ + δ + L_tilde_derivative_term - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2
        # and we are decomposing this into the terms involving (I with L_1) AND (L_2) -- otherwise it will perform elementwise operation 
        # on two banded matrices which takes extreme amount of memory that are not needed.
        residual[1:M] = (ρ + δ + L_tilde_derivative_term - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:M] # system of ODEs (eq:28)
        residual[1:M] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:M] 
        residual[1:M] .-= (υ^2/2)*L_2*u[1:M] 
        residual[1:M] .-= du[1:M]
        residual[1:M] .-= π_tilde
        residual[M+1] = u[1] + x - dot(ω, u[1:M]) # residual (eq:25)
        residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31)
    end

    u0 = [p.v_T; p.g_T; p.z_hat_T]
    du0 = zeros(M+2)

    return (dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), differential_vars = [trues(M); false; false], p), stationary_equilibrium = stationary_equilibrium, f! = f!, S = get_S, L_tilde = get_L_tilde)
end