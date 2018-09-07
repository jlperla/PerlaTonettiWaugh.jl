function solve_dynamics(params, settings, d_0, d_T)
    @unpack δ = params
    @unpack z = settings 
    M = length(z)

    # Compute the stationary solution at t = 0 and t = T
    params_0 = merge(params, (d = d_0,)) # parameters to be used at t = 0
    params_T = merge(params, (d = d_T,)) # parameters to be used at t = T

    stationary_sol_0 = stationary_numerical(params_0, z) # solution at t = 0
    Ω_0 = stationary_sol_0.Ω

    stationary_sol_T = stationary_numerical(params_T, z) # solution at t = T
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    Ω_T = stationary_sol_T.Ω

    # compute the resulting end time and function of Ω
    T = (log(Ω_0) - log(Ω_T)) / δ
    Ω(t) = t < T ? Ω_0 * exp(-δ*t) : Ω_T
    E(t) = t >= T ? δ : 0
    E(t) = 1 # TODO: remove this later to see if discontinuity is resolved

    # define the corresponding DAE problem
    p = get_p(params_T, stationary_sol_T, settings, Ω, T)
    dae = PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)

    # solve solutions
    sol = DifferentialEquations.solve(dae.dae_prob, callback = dae.callback) # solve!
    residuals = calculate_residuals_dae(dae.dae_prob.f, deepcopy(sol.du), deepcopy(sol.u), p, sol.t)

    v = map(u -> u[1:M], sol.u)
    g = map(u -> u[M+1], sol.u)
    z_hat = map(u -> u[M+2], sol.u)

    return (v = v, g = g, z_hat = z_hat, t = sol.t, 
            p = p, sol = sol, residuals = residuals)
end

# Implementation of the full model with time-varying objects, represented by DAE
function PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)
    # Unpack params and settings. 
    @unpack z = settings 
    M = length(z)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, θ, δ, χ, N, ζ, ρ = p 

    function stationary_equilibrium(g, z_hat, E, Ω, t)
        S = (g - μ - θ * υ^2/2)
        L_tilde = Ω * ((N-1) * z_hat^(-θ)*κ + ζ*θ*S + ζ*E*δ / χ)
        z_bar = Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        return (S = S, L_tilde = L_tilde, z_bar = z_bar, π_min = π_min)
    end

    callback = SavingCallback((u,t,integrator)->(t, stationary_equilibrium(u[M+1], u[M+2], E(t), Ω(t), t).L_tilde), 
                                p.saved_values, 
                                tdir = -1) # need to compute D_t L(t)

    # Dynamic calculations, defined for each time ∈ t.
    function f!(residual,du,u,p,t)
        residual[:] = zeros(M+2)

        # Carry out calculations.
        v = u[1:M]
        g = u[M+1]
        z_hat = u[M+2]
        @unpack S, L_tilde, z_bar, π_min = stationary_equilibrium(g, z_hat, E(t), Ω(t), t)
        π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= log(z_hat))) - (N-1)*κ*exp(-(σ-1)*z)*(z >= log(z_hat))
        # π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= z_hat)) - (N-1)*κ*exp(-(σ-1)*z)*(z >= z_hat)
        x = ζ
        # compute the derivative of L_tilde
        values_future = p.saved_values.saveval
        L_tilde_derivative = 0 # default value
        forward_index = findlast(x -> x[1] > t, values_future)
        if (forward_index != nothing;)
            if (forward_index > 0)
            t_forward = values_future[forward_index][1]
            L_tilde_t_forward = values_future[forward_index][2]
            L_tilde_t_derivative = (L_tilde_t_forward - L_tilde) / (t_forward - t)
            end
        end

        # form the DAE at t
        ρ_tilde = ρ + δ + L_tilde_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2)
        A_t = ρ_tilde*I - (μ - g + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2        
        residual[1:M] = A_t * v - π_tilde.(z) # system of ODEs (eq:28)
        residual[M+1] = v[1] + x - dot(ω, v) # residual (eq:25)
        residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31) 
        residual[1:M] .-= du[1:M]    
    end

    u0 = [p.v_T; p.g_T; p.z_hat_T]
    du0 = zeros(M+2)

    return (dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), differential_vars = [trues(M); false; false], p),
    callback = callback)
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

# calculate the residual at each time point; each row represents t in ts
function calculate_residuals_dae(f!, du, u, p, ts)    
    residuals = zeros(length(ts), length(u[1]))
    
    for (i, t) in enumerate(ts)
        residual = zeros(length(u[1]))
        f!(residual, du[i], u[i], p, t)
        residuals[i,:] = residual 
    end

    return residuals
end
