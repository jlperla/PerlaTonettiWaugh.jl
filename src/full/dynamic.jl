# Kernel function for main method.
function f!(residual,du,u,p,t)
    @unpack results, parameters, settings, L_1, L_2, ω, E, Ω, Ξ₁ = p 
    @unpack z, T = settings 
    @unpack ρ, δ, σ, μ, υ, κ, d, ζ = parameters
    residual .= 0
    P = length(residual) - 2 # note u[1:P]  = v(t) in the solution iterators
    g = u[P+1]
    z_hat = u[P+2]
    x = ζ
    Ω_t = Ω(t)
    E_t = E(t)

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S_t, L_tilde_t, z_bar, π_min, π = static_equilibrium(Ξ₁, u[1], g, z_hat, E_t, Ω_t, z, parameters)
    L_tilde_log_derivative = 0.0 # Default to Float literal.
    if (t < T)
        t_forward = results[:t][end]
        L_tilde_forward = results[:L_tilde][end]
        L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde_t))/(t_forward - t) # See note under (40)
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
function solve_dynamics(parameters, sol, settings, Ω, E; detailed_solution = true)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z, z_ex, T, tstops = settings
    @assert γ ≈ 1 # These are the only supported parameters for the transition dynamics at this point
    @assert η == 0
    T >= 50 ? nothing : @warn "T < 50 can give erratic behavior."

    # unpack the stationary solution 
    v_T, g_T, z_hat_T, L_tilde_T, Ω_T = (sol.v_tilde, sol.g, sol.z_hat, sol.L_tilde, sol.Ω)

    # Define the results data frame we'll be using and push the stationary onto it.
    results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_1 = v_T[1], L_tilde = L_tilde_T)

    # Define intermediate quantitities.
    P = length(z)
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)

    # initial conditions and parameters for the solver
    u0 = [v_T; g_T; z_hat_T]
    du0 = zeros(P+2)

    p = (results = results, parameters = parameters, settings = settings, 
        L_1 = L_1, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, E = E, Ω = Ω)

    dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false])

    # called at end of each time-step in the DAE solver (primiarly to calculate L_tilde for derivatives)
    function cb_aux(u, t, integrator)
        g_t = u[P+1]
        z_hat_t = u[P+2]
        Ω_t = Ω(t)
        E_t = E(t)
        v_1_t = u[1]
        S_t = S(g_t, parameters)
        L_tilde_t = L_tilde(g_t, z_hat_t, Ω_t, E_t, S_t, parameters)
        push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_1 = v_1_t, L_tilde = L_tilde_t))
    end

    # Solve that DAE passing in the callbacks to store "future" L_tilde
    cb = FunctionCallingCallback(cb_aux, tdir = -1, func_start = false) # Callback object.
    if (tstops == nothing)
        sol = DifferentialEquations.solve(dae_prob, callback = cb)
    else
        sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops)
    end

    # postprocess results and return 
    sort!(results) # sort results to be forward in time 
    # define auxiliary functions 
    entry_residual(v_1) = Ξ₁*v_1 - ζ*(1-χ)/χ # (56)  
    z_bar(Ω_t, z_hat) = ((Ω_t * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)))^(1/(σ-1))) # (37)
    π_min(L_tilde_t, z_bar) = (1 - L_tilde_t) / ((σ-1)*z_bar^(σ-1)) # (38)
    # define interpolates
    g_interpolated(t) = (sol(t))[P+1]
    z_hat_interpolated(t) = (sol(t))[end]
    L_tilde_interpolated(t) = L_tilde(g_interpolated(t), z_hat_interpolated(t), Ω(t), E(t), S(g_interpolated(t), parameters), parameters)
    z_bar_interpolated(t) = z_bar(Ω(t), z_hat_interpolated(t))
    c_interpolated(t) = c(L_tilde_interpolated(t), Ω(t), z_bar_interpolated(t))
    π_min_interpolated(t) = π_min(L_tilde_interpolated(t), z_bar_interpolated(t))
    # use interpolates to compute quantities 
    log_c(t) = log(c_interpolated(t))
    log_c_T = log(c(L_tilde_T, Ω_T, z_bar(Ω_T, z_hat_T))) # this is a scalar 
    log_M(t) = quadgk(g_interpolated, 0, t)[1]
    U(t) = quadgk(τ -> exp(-ρ*τ)*(log_M(t+τ) + log_c(t+τ)), 0, (T-t))[1] + exp(-ρ*(T-t))*(g_T + ρ*(log_c_T + g_T * T))/(ρ^2) # (C.83)
    π_rat(t) = θ*(1-z_hat_interpolated(t)^(-θ+σ-1))/(θ-σ+1) + (1+(N-1)*d^(1-σ))*θ*z_hat_interpolated(t)^(-θ+σ-1)/(θ-σ+1)+(N-1)*κ*z_hat_interpolated(t)^(-θ)/π_min_interpolated(t) # (61)

    # Add these quantities to the DataFrame.
    results = @transform(results, entry_residual = entry_residual.(:v_1))

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
        results = @transform(results, entry_residual = entry_residual.(:v_1))

        results = @transform(results, λ_ii = λ_ii.(:z_hat, Ref(parameters))) # λ_ii column.
        results = @transform(results, S = S.(:g, Ref(parameters))) # S column.
        results = @transform(results, z_bar = z_bar.(:Ω, :z_hat)) # z_bar column.
        results = @transform(results, c = c.(:L_tilde, :Ω, :z_bar)) # c column.
        results = @transform(results, π_min = π_min.(:L_tilde, :z_bar)) # π_min column.
        results = @transform(results, log_M = log_M.(:t)) # log_M column
        results = @transform(results, U = U.(:t)) # U column
        results = @transform(results, π_rat = π_rat.(:t)) # π_rat column
        results = @transform(results, L_tilde_a = L_tilde_a.(:Ω, :S, Ref(parameters)))
        results = @transform(results, L_tilde_x = L_tilde_x.(:z_hat, :Ω, Ref(parameters)))
        results = @transform(results, L_tilde_E = L_tilde_E.(:Ω, :E, Ref(parameters)))
        results = @transform(results, w = w.(:z_bar, Ref(parameters)))

        # logic for r
        results.r = ones(Float64, nrow(results)) # filler, to be overwritten
        for i in 1:nrow(results)
            t = results[:t][i]
            c = results[:c][i]
            g = results[:g][i]
            log_c_forward = (i < nrow(results)) ? (log(results[:c][i+1]) - log(c))/(results[:t][i+1] - t) : 0.0
            if (i != nrow(results)) 
                @assert results[:t][i+1] > t # ensure that differencing is actually forward
            end
            results.r[i] = ρ + γ*(g + log_c_forward) # (C.56)
        end

    end

    # The results, raw DAE solution, and DAE problem (f!, static_equilibrium, etc.) objects.
    return (results = results, sol = sol, p = p,
            U = U, c = c, Ω = Ω, log_M = log_M, log_c = log_c)
end
