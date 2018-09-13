function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)
    @unpack δ, N, σ, θ, d = params_T
    @unpack z, tstops, Δ_E = settings 
    M = length(z)

    # define E(t) based on FD
    E(t) = t < (T - Δ_E) ? (log(Ω(t+Δ_E)) - log(Ω(t)))/Δ_E + δ : δ
    
    # define the corresponding DAE problem
    p = get_p(params_T, stationary_sol_T, settings, Ω, T)
    dae = PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)

    # solve solutions
    sol = DifferentialEquations.solve(dae.dae_prob, callback = dae.callback, tstops = tstops) # solve!
    @unpack u, du, t = sol

    residuals = zeros(length(t), length(u[1]))
    equilibriums = []
    for (i, t) in enumerate(t)
        # compute residual at t
        residual = zeros(length(u[1]))
        dae.dae_prob.f(residual, du[i], u[i], p, t)
        residuals[i,:] = residual 
        # compute stationary equilibrium at t 
        v_1 = u[i][1]
        g = u[i][M+1]
        z_hat = u[i][M+2]
        # compute λ_ii and c
        equilibrium = dae.stationary_equilibrium(v_1, g, z_hat, E(t), Ω(t), t)
        λ_ii = 1 / (1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ))
        c = (θ / (1-σ+θ))^(1/(σ-1))*(1-equilibrium.L_tilde)*Ω(t)^(1/(σ-1))*λ_ii^(1/(1-σ))
        # TODO: add log_M and U later
        equilibrium = merge(equilibrium, (λ_ii = λ_ii, c = c,))
        push!(equilibriums, equilibrium)
    end

    # extract solutions
    v = map(u -> u[1:M], sol.u)
    g = map(u -> u[M+1], sol.u)
    z_hat = map(u -> u[M+2], sol.u)
    S = map(eq -> eq.S, equilibriums)
    L_tilde = map(eq -> eq.L_tilde, equilibriums)
    z_bar = map(eq -> eq.z_bar, equilibriums)
    π_min = map(eq -> eq.π_min, equilibriums)
    π_tilde = map(eq -> eq.π_tilde, equilibriums)
    λ_ii = map(eq -> eq.λ_ii, equilibriums)
    c = map(eq -> eq.c, equilibriums)
    entry_residual = map(eq -> eq.entry_residual, equilibriums)

    return (v = v, g = g, z_hat = z_hat, 
            S = S, L_tilde = L_tilde, z_bar = z_bar, π_min = π_min, π_tilde = π_tilde,
            λ_ii = λ_ii, c = c, entry_residual = entry_residual,
            t = t, p = p, sol = sol, residuals = residuals, equilibriums = equilibriums)
end

# Implementation of the full model with time-varying objects, represented by DAE
function PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)
    # Unpack params and settings. 
    @unpack z = settings 
    M = length(z)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, θ, δ, χ, N, ζ, ρ = p 

    function stationary_equilibrium(v_1, g, z_hat, E, Ω, t)
        S = θ * (g - μ - θ * υ^2/2)
        L_tilde = Ω * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E / χ))
        z_bar = Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= log(z_hat))) - (N-1)*κ*exp(-(σ-1)*z)*(z >= log(z_hat))
        # π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= z_hat)) - (N-1)*κ*exp(-(σ-1)*z)*(z >= z_hat)
        π_tilde = π_tilde.(z)
        entry_residual = v_1 - ζ * (1-χ) / χ
        return (S = S, L_tilde = L_tilde, z_bar = z_bar, 
                π_min = π_min, π_tilde = π_tilde,
                entry_residual = entry_residual)
    end

    callback = SavingCallback((u,t,integrator)->(t, stationary_equilibrium(u[1], u[M+1], u[M+2], E(t), Ω(t), t).L_tilde), 
                                p.saved_values, 
                                tdir = -1) # need to compute D_t L(t)

    # Dynamic calculations, defined for each time ∈ t.
    function f!(residual,du,u,p,t)
        residual[:] = zeros(M+2)

        # Carry out calculations.
        v = u[1:M]
        g = u[M+1]
        z_hat = u[M+2]
        @unpack S, L_tilde, z_bar, π_min, π_tilde = stationary_equilibrium(v[1], g, z_hat, E(t), Ω(t), t)
        x = ζ
        # compute the derivative of L_tilde
        values_future = p.saved_values.saveval
        L_tilde_derivative_term = (params_T.γ - 1) * g # default value
        forward_index = findlast(x -> x[1] > t, values_future)
        if ((T - t) > 1e-3 && forward_index != nothing;) # use callbacks only if t is well-separated from T 
            if (forward_index > 0)
            t_forward = values_future[forward_index][1]
            L_tilde_t_forward = values_future[forward_index][2]
            L_tilde_derivative_term = (log(1 - L_tilde_t_forward) - log(1 - L_tilde)) / (t_forward - t)
            end
        end

        # form the DAE at t
        ρ_tilde = ρ + δ + L_tilde_derivative_term - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2)
        A_t = ρ_tilde*I - (μ - g + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2        
        residual[1:M] = A_t * v - π_tilde # system of ODEs (eq:28)
        residual[M+1] = v[1] + x - dot(ω, v) # residual (eq:25)
        residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31) 
        residual[1:M] .-= du[1:M]    
    end

    u0 = [p.v_T; p.g_T; p.z_hat_T]
    du0 = zeros(M+2)

    return (dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), differential_vars = [trues(M); false; false], p),
            callback = callback, stationary_equilibrium = stationary_equilibrium)
end

# return the parameters and functions needed to define dynamics
function get_p(params_T, stationary_sol_T, settings, Ω, T)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T
    @unpack z = settings 
    M = length(z)

    # Unpack the stationary solution.
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    Ω_T = stationary_sol_T.Ω

    # Quadrature weighting
    ω = ω_weights(z, θ, σ-1)  # TODO: compute integral for eq 30

    # Discretize the operator. 
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # L_1_minus ≡ L_1 is the only one we use. 

    # Bundle as before.
    p = (L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ, 
        ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, ω = ω, Ω = Ω,
        v_T = v_T, g_T = g_T, z_hat_T = z_hat_T, Ω_T = Ω_T,
        saved_values = SavedValues(Float64, Tuple{Float64,Float64})) #Named tuple for parameters.
    return p
end
