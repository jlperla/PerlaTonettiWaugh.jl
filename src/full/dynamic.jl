function solve_dynamic_full(params, settings, d_0, d_T)
    @unpack δ = params
    @unpack z = settings 
    M = length(z)

    # Compute the stationary solution at t = 0 and t = T
    params_0 = merge(params, @NT(d = d_0)) # parameters to be used at t = 0
    params_T = merge(params, @NT(d = d_T)) # parameters to be used at t = T

    stationary_sol_0 = stationary_numerical(params_0, z) # solution at t = 0
    Ω_0 = stationary_sol_0.Ω

    stationary_sol_T = stationary_numerical(params_T, z) # solution at t = T
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    Ω_T = stationary_sol_T.Ω

    # compute the resulting end time and function of Ω
    T = (log(Ω_0) - log(Ω_T)) / δ
    Ω = Ω(t) = t < T ? Ω_0 * exp(-δ*t) : Ω_T

    # define the corresponding DAE problem
    p = get_p(params_T, stationary_sol_T, settings, Ω, T)
    dae = fullDAE(params_T, stationary_sol_T, settings, Ω, T, p)

    # solve solutions
    @time sol = DifferentialEquations.solve(dae.dae_prob, callback = dae.callback) # solve!
    @time residuals = calculate_residuals(sol.du, sol.u, p, sol.t)

    return @NT(sol = sol, p = p, residuals = residuals)
end

function calculate_residuals(du, u, p, ts)
    # Calculate the residual at each time point
    residuals = zeros(length(ts), length(u[1]))
    
    for (i, t) in enumerate(ts)
        residuals[i,:] = calculate_residual_t(du[i], u[i], p, t)
    end

    return residuals
end


# Implementation of the full model with time-varying objects, represented by DAE
function fullDAE(params_T, stationary_sol_T, settings, Ω, T, p)
    # Unpack params and settings. 
    @unpack z = settings 
    M = length(z)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, θ, δ, χ, N, ζ, ρ = p 

    function stationary_equilibrium(g, z_hat, Ω, t)
        S = (g - μ - θ * υ^2/2)
        E = t >= T ? δ : 0
        E = 1 # TODO: remove this later to see if discontinuity is resolved
        L_tilde = Ω * ((N-1) * z_hat^(-θ)*κ + ζ*θ*S + ζ*E*δ / χ)
        z_bar = Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        return @NT(S = S, E = E, L_tilde = L_tilde, z_bar = z_bar, π_min = π_min)
    end

    callback = SavingCallback((u,t,integrator)->(t, stationary_equilibrium(u[M+1], u[M+2], Ω(t), t).L_tilde), 
                                p.saved_values, 
                                tdir = -1) # need to compute D_t L(t)
    
    # Dynamic calculations, defined for each time ∈ t.  
    function f!(resid,du,u,p,t) 
        resid[:] = zeros(u)
        
        # Carry out calculations. 
        v = u[1:M]
        g = u[M+1]
        z_hat = u[M+2]
        @unpack S, E, L_tilde, z_bar, π_min = stationary_equilibrium(g, z_hat, Ω(t), t)
        π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= log(z_hat))) - (N-1)*κ*exp(-(σ-1)*z)*(z >= log(z_hat))
        # π_tilde(z) = π_min * (1+(N-1)*d^(1-σ)*(z >= z_hat)) - (N-1)*κ*exp(-(σ-1)*z)*(z >= z_hat)
        x = ζ
        # compute the derivative of L_tilde
        values_future = p.saved_values.saveval
        L_tilde_derivative = 0 # default value
        forward_index = findlast(x -> x[1] > t, values_future)
        if (forward_index > 0)
            t_forward = values_future[forward_index][1]
            L_tilde_t_forward = values_future[forward_index][2]
            L_tilde_t_derivative = (L_tilde_t_forward - L_tilde) / (t_forward - t)
        end

        # form the DAE at t
        ρ_tilde = ρ + δ + L_tilde_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2)
        A_t = ρ_tilde*I - (μ - g + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2        
        resid[1:M] = A_t * v - π_tilde.(z) # system of ODEs (eq:28)
        resid[M+1] = v[1] + x - dot(ω, v) # residual (eq:25)
        resid[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31) 
        resid[1:M] .-= du[1:M]    
    end

    u = [p.v_T; p.g_T; p.z_hat_T]
    du = zeros(M+2)
    resid_M2 = zeros(M+2)

    return @NT(dae_prob = DAEProblem(f!, resid_M2, u, (T, 0.0), differential_vars = [trues(M); false; false], p),
    callback = callback)
end

function stationary_equilibrium(g, z_hat, Ω_static, t)

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
    p = @NT(L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ, 
        ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, ω = ω, Ω = Ω,
        v_T = v_T, g_T = g_T, z_hat_T = z_hat_T, Ω_T = Ω_T,
        saved_values = SavedValues(Float64, Tuple{Float64,Float64})) #Named tuple for parameters.
    return p
end

function calculate_residual_t(du, u, p, t)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, Ω = p 
    resid = zeros(u)
    
    # Carry out calculations. 
    v_t = u[1:M]
    g_t = u[M+1]
    z_hat_t = u[M+2]

    @unpack x_t, π_min_t, π_tilde_t_by_z, ρ_tilde_t = get_static_vals(p, t, v_t, g_t, z_hat_t)

    A_t = ρ_tilde_t*I - (μ - g_t + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2        
    resid[1:M] = A_t * v_t - π_tilde_t_by_z # system of ODEs (eq:28)
    resid[M+1] = v_t[1] + x_t - dot(ω, v_t) # residual (eq:25)
    resid[M+2] = z_hat_t^(σ-1) - κ * d^(σ-1) / π_min_t # export threshold (eq:31) 
    resid[1:M] .-= du[1:M]

    return resid
end

function get_L_tilde_t(p, t, g_t, z_hat_t)
    @unpack N, T, θ, κ, ζ, δ, χ, Ω, μ, υ = p
    is_t_over_T = ifelse(t >= T, 1, 0)
    is_t_over_T = 1 # FIXIT: later see if removing this line still makes no discontinuity
    return Ω(t) * ((N-1) * z_hat_t^(-θ)*κ + ζ*θ*(g_t - μ - θ * υ^2/2) + is_t_over_T * ζ * δ / χ)
end

function get_static_vals(p, t, v_t, g_t, z_hat_t)
    @unpack z, N, T, θ, σ, κ, ζ, d, ρ, δ, υ, μ, χ, Ω, saved_values = p

    is_z_over_log_z_hat_t = z_val -> ifelse(z_val >= log(z_hat_t), 1, 0)

    # compute L_tilde_t and its derivative
    L_tilde_t = get_L_tilde_t(p, t, g_t, z_hat_t)
    values_future = saved_values.saveval
    L_tilde_t_derivative = 0 # default value
    forward_index = findlast(x -> x[1] > t, values_future)
    if (forward_index > 0)
        t_forward = values_future[forward_index][1]
        L_tilde_t_forward = values_future[forward_index][2]
        L_tilde_t_derivative = (L_tilde_t_forward - L_tilde_t) / (t_forward - t)
    end

    # compute other variables evaluated at t (see subsection 2.5 in the note)
    x_t = ζ
    z_bar_t_term = Ω(t) * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat_t^(σ-1-θ))
    π_min_t = (1 - L_tilde_t) / ((σ-1)*z_bar_t_term)
    π_tilde_t_map = z_val -> (π_min_t * (1+(N-1)*d^(1-σ)*is_z_over_log_z_hat_t(z_val)) - (N-1)*κ*exp(-(σ-1)*z_val)*is_z_over_log_z_hat_t(z_val))

    π_tilde_t_by_z = map(π_tilde_t_map, z)
    r_tilde_t = ρ + δ + L_tilde_t_derivative
    ρ_tilde_t = r_tilde_t - (σ-1) * (μ - g_t + (σ-1)*υ^2 / 2)

    return @NT(x_t = x_t, π_min_t = π_min_t, π_tilde_t_by_z = π_tilde_t_by_z, ρ_tilde_t = ρ_tilde_t)
end
