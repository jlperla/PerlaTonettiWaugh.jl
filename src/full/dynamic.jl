T_Δ_MIN = 1e-03

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
    Ω = get_Ω(Ω_T, δ, T)

    # define the corresponding DAE problem
    p = get_p(params_T, stationary_sol_T, settings, Ω, T)
    dae_prob = fullDAE(params_T, stationary_sol_T, settings, Ω, T, p)

    # solve solutions
    tstops = 0:T_Δ_MIN:T # ensure that time grids are fine enough

    callback = SavingCallback((u,t,integrator)->(t, get_L_tilde_t(p, t, u[M+1], p.map_z_hat_t(u[M+2]))), 
                p.saved_values, 
                tdir = -1) # need to compute D_t L(t)
    sol = DifferentialEquations.solve(dae_prob, callback = callback, tstops = tstops) # solve!

    return @NT(sol = sol, p = p)
end

# Implementation of the full model with time-varying objects, represented by DAE
function fullDAE(params_T, stationary_sol_T, settings, Ω, T, p)
    # Unpack params and settings. 
    @unpack z = settings 
    M = length(z)
    
    # Dynamic calculations, defined for each time ∈ t.  
    function f!(resid,du,u,p,t)
        resid[:] = calculate_residual_t(du, u, p, t)
    end

    u = [p.v_T; p.g_T; p.map_z_hat_t_inverse(p.z_hat_T)]
    du = zeros(M+2)
    resid_M2 = zeros(M+2)

    return DAEProblem(f!, resid_M2, u, (T, 0.0), differential_vars = [trues(M); false; false], p)
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
    map_z_hat_t = bridge(ℝ, Segment(1, 100))
    map_z_hat_t_inverse = val -> inverse(map_z_hat_t, val) 
    p = @NT(L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ, 
        ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, ω = ω, Ω = Ω,
        v_T = v_T, g_T = g_T, z_hat_T = z_hat_T, Ω_T = Ω_T,
        map_z_hat_t = map_z_hat_t, 
        map_z_hat_t_inverse= map_z_hat_t_inverse,
        saved_values = SavedValues(Float64, Tuple{Float64,Float64})) #Named tuple for parameters.
    return p
end

function calculate_residual_t(du, u, p, t)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, Ω, map_z_hat_t = p 
    resid = zeros(u)
    
    # Carry out calculations. 
    v_t = u[1:M]
    g_t = u[M+1]
    z_hat_t = map_z_hat_t(u[M+2])

    @unpack x_t, π_min_t, π_tilde_t_by_z, ρ_tilde_t = get_static_vals(p, t, v_t, g_t, z_hat_t)

    A_t = ρ_tilde_t*I - (μ - g_t + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2        
    resid[1:M] = A_t * v_t - π_tilde_t_by_z # system of ODEs (eq:28)
    resid[M+1] = v_t[1] + x_t - dot(ω, v_t) # residual (eq:25)
    resid[M+2] = z_hat_t^(σ-1) - κ * d^(σ-1) / π_min_t # export threshold (eq:31) 
    resid[1:M] .-= du[1:M]

    return resid
end

function get_Ω(Ω_T, δ, T)
    Ω_0 = Ω_T * exp(δ*T)
    function Ω(t)
        if (t < T)
            return Ω_0 * exp(-δ*t)
        end
        return Ω_T
    end
    return Ω
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
    # if (forward_index > 0)
    if (T - t > T_Δ_MIN)
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
