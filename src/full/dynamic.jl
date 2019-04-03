# Kernel function for main method.
function f!(residual,du,u,p,t)
    @unpack ζ, Ω, E, static_equilibrium, T, results, ρ, δ, σ, μ, υ, L_1, L_2, ω, κ, d, Ξ₁ = p
    residual .= 0
    P = length(residual) - 2 # note u[1:P]  = v(t) in the solution iterators
    g = u[P+1]
    z_hat = u[P+2]
    x = ζ
    Ω_t = Ω(t)
    E_t = E(t)

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S_t, L_tilde_t, z_bar, π_min, π = static_equilibrium(u[1], g, z_hat, E_t, Ω_t)
    L_tilde_log_derivative = 0.0 # Default to Float literal.
    if (t < T)
        t_forward = results[:t][end]
        L_tilde_forward = results[:L_tilde][end]
        L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde_t))/(t_forward - t) # See note under (34)
    end

    # Combine (40) with (52) to get A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2
    # Then building the A(t) * v(t) - v'(t) residual directly for the non-algebraic equations
    residual[1:P] = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:P] # (40) into (52) then (53)
    residual[1:P] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:P] # (52) 2nd term in (53)
    residual[1:P] .-= (υ^2/2)*L_2*u[1:P] # (52) third term in (53)
    residual[1:P] .-= π # (53) final term
    residual[1:P] .-= du[1:P] # (53) subtracting the v'(t) to form the residual
    residual[P+1] = Ξ₁*u[1] + ζ - dot(ω, u[1:P])  # (54)
    residual[P+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min  # (55)
end

# Calculate the transition dynamics given a fixed Ω(t), E(t) function
# To calculate the full equilibrium, you can solve until Ω(t), E(t) fulfill free-entry conditions
function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; detailed_solution = true)

    if (T < settings.T_U_bar)
        throw("Terminal time `T` should be large enough so that T >= settings.T_U_bar is satisfied.")
    end

    # Unpack arguments
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T
    @unpack z, z_ex, T_U_bar, tstops = settings
    @assert γ ≈ 1 # These are the only supported parameters for the transition dynamics at this point
    @assert η == 0

    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    L_tilde_T = stationary_sol_T.L_tilde
    Ω_T = stationary_sol_T.Ω

    # Define the results data frame we'll be using and push the stationary onto it.
    results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_1 = v_T[1], L_tilde = L_tilde_T)

    # Define intermediate quantitities.
    P = length(z)
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂(z_ex, bc)
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)

    # Define the auxiliary functions for the DAE problem.
    S(g) = θ * (g - μ - θ * υ^2/2)  # (32)
    L_tilde(S, z_hat, E_t, Ω_t) = Ω_t * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E_t / χ))  # (33)

    function static_equilibrium(v_1, g, z_hat, E_t, Ω_t)
        S_t = S(g)
        L_tilde_t = L_tilde(S_t, z_hat, E_t, Ω_t)
        z_bar = (Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)))^(1/(σ-1))  # (37)
        w = ((σ-1)/σ)*z_bar # (C.13)
        π_min = (1 - L_tilde_t) / ((σ-1)*z_bar^(σ-1))  # (38)

        i_vectorized = z .>= log(z_hat) # Vectorized indicator function
        π = π_min * (1.0.+(N-1)*d^(1-σ)*i_vectorized) - (N-1)*κ*exp.(-(σ-1).*z).*i_vectorized  # (39)
        entry_residual = Ξ₁*v_1 - ζ * (1-χ) / χ  # (56)
        return (S_t = S_t, L_tilde_t = L_tilde_t, z_bar = z_bar, π_min = π_min, π = π, entry_residual = entry_residual,
                w = w)
    end

    # initial conditions and parameters for the solver
    u0 = [v_T; g_T; z_hat_T]
    du0 = zeros(P+2)

    p = (ζ = ζ, Ω = Ω, E = E, static_equilibrium = static_equilibrium, T = T,
        results = results, ρ = ρ, δ = δ, σ = σ, μ = μ, υ = υ, L_1 = L_1, L_2 = L_2,
        ω = ω, κ = κ, d = d, Ξ₁ = Ξ₁)
    dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false])

    # called at end of each time-step in the DAE solver (primiarly to calculate L_tilde for derivatives)
    function cb_aux(u, t, integrator)
        g_t = u[P+1]
        z_hat_t = u[P+2]
        Ω_t = Ω(t)
        E_t = E(t)
        v_1_t = u[1]
        S_t = S(g_t)
        L_tilde_t = L_tilde(S_t, z_hat_t, E_t, Ω_t)
        push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_1 = v_1_t, L_tilde = L_tilde_t))
    end


    # Solve that DAE passing in the callbacks to store "future" L_tilde
    cb = FunctionCallingCallback(cb_aux, tdir = -1, func_start = false) # Callback object.
    if (tstops == nothing)
        sol = DifferentialEquations.solve(dae_prob, callback = cb)
    else
        sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops)
    end

    # add calculations to the DataFrame
    results = sort!(results) # sort results to be forward in time

    gen_λ_ii(z_hat) = 1 / (1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (57)
    gen_c(L_tilde, Ω, z_bar, S) = (1 - L_tilde)*z_bar # (58)
    gen_S = S
    gen_z_bar(Ω_t, z_hat) = ((Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)))^(1/(σ-1)))^(1/(σ-1)) # (37)
    gen_π_min(L_tilde_t, z_bar) = (1 - L_tilde_t) / ((σ-1)*z_bar^(σ-1)) # (38)
    gen_entry_residual(v_1) = Ξ₁*v_1 - ζ*(1-χ)/χ # (56)
    gen_L_tilde_adopt(Ω, S) = Ω * ζ * S # (36)
    gen_L_tilde_export(Ω, z_hat) = Ω * ((N-1)*z_hat^(-θ))*κ # (34)
    gen_L_tilde_entrycost(Ω, E) = Ω * ζ * E / χ # (35)
    gen_w(z_bar) = ((σ-1)/σ) * z_bar # (C.13)

    # Add these quantities to the DataFrame.
    results = @transform(results, entry_residual = gen_entry_residual.(:v_1)) # entry_residual column
    log_c_T = log(gen_c(L_tilde_T, Ω_T, gen_z_bar(Ω_T, z_hat_T), S(g_T)))

    g_interpolated(t) = (sol(t))[P+1]
    z_hat_interpolated(t) = (sol(t))[end]
    L_tilde_interpolated(t) = L_tilde(S(g_interpolated(t)), z_hat_interpolated(t), E(t), Ω(t))
    z_bar(t) = gen_z_bar(Ω(t), z_hat_interpolated(t))
    U(t) = quadgk(τ -> exp(-ρ*τ)*(log_M(t+τ) + log_c(t+τ)), 0, (T-t))[1] + exp(-ρ*(T-t))*(g_T + ρ*(log_c_T + g_T * T))/(ρ^2)  # (C.83)
    c(t) = gen_c(L_tilde_interpolated(t), Ω(t), z_bar(t), S(t))
    log_M(t) = quadgk(g_interpolated, 0, t)[1]
    log_c(t) = log(gen_c(L_tilde_interpolated(t), Ω(t), gen_z_bar(Ω(t), z_hat_interpolated(t)), S(g_interpolated(t))))
    π_min(t) = gen_π_min(L_tilde_interpolated(t), z_bar(t))
    π_rat(t) = θ*(1-z_hat_interpolated(t)^(-θ+σ-1))/(θ-σ+1) + (1+(N-1)*d^(1-σ))*θ*z_hat_interpolated(t)^(-θ+σ-1)/(θ-σ+1)+(N-1)*κ*z_hat_interpolated(t)^(-θ)/π_min(t) # (61)

    # These extra calculations are best skipped if running the optimizer
    if (detailed_solution)
        ts = (T >= 20) ? sort(unique([collect(0.0:0.5:20); collect(20.0:T); results.t])) : results.t # add more time steps for the plotting/etc.
        results = DataFrame(t = ts)

        results = @transform(results, g = g_interpolated.(:t))
        results = @transform(results, z_hat = z_hat_interpolated.(:t))
        results = @transform(results, Ω = Ω.(:t))
        results = @transform(results, E = E.(:t))
        results = @transform(results, v_1 = (t -> (sol(t))[1]).(:t))
        results = @transform(results, L_tilde = L_tilde_interpolated.(:t))
        results = @transform(results, entry_residual = gen_entry_residual.(:v_1))

        results = @transform(results, λ_ii = gen_λ_ii.(:z_hat)) # λ_ii column.
        results = @transform(results, S = gen_S.(:g)) # S column.
        results = @transform(results, z_bar = gen_z_bar.(:Ω, :z_hat)) # z_bar column.
        results = @transform(results, c = gen_c.(:L_tilde, :Ω, :z_bar, :S)) # c column.
        results = @transform(results, π_min = gen_π_min.(:L_tilde, :z_bar)) # π_min column.
        results = @transform(results, log_M = log_M.(:t)) # log_M column
        results = @transform(results, U = U.(:t)) # U column
        results = @transform(results, π_rat = π_rat.(:t)) # π_rat column
        results = @transform(results, L_tilde_a = gen_L_tilde_adopt.(:Ω, :S))
        results = @transform(results, L_tilde_x = gen_L_tilde_export.(:Ω, :z_hat))
        results = @transform(results, L_tilde_E = gen_L_tilde_entrycost.(:Ω, :E))
        results = @transform(results, w = gen_w.(:z_bar))

        # logic for r
        results.r = ones(Float64, nrow(results)) # filler, to be overwritten
        for i in 1:nrow(results)
            t = results[:t][i]
            c = results[:c][i]
            g = results[:g][i]
            log_c_forward = (i < nrow(results)) ? (log(results[:c][i+1]) - log(c))/(results[:t][i+1] - t) : 0.0
            i == nrow(results) || @assert results[:t][i+1] > t # ensure that differencing is actually forward
            results.r[i] = ρ + γ*(g + log_c_forward) # (C.56)
        end

    end

    # The results, raw DAE solution, and DAE problem (f!, static_equilibrium, etc.) objects.
    return (results = results, sol = sol, p = p, static_equilibrium = static_equilibrium,
          U = U, c = c, Ω = Ω, log_M = log_M, log_c = log_c)
end
