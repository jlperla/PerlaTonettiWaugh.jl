
function solve_dynamic_full(params, settings)
    @unpack z, T = settings 
    M = length(z)
    tstops = 0:1e-03:T # TODO: might need finer grids
    dae_prob = fullDAE(params, settings)
    
    DifferentialEquations.solve(dae_prob, tstops=tstops)
end

# Implementation of the full model with time-varying objects, represented by DAE
function fullDAE(params, settings)
    # Unpack params and settings. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack z, T = settings 
    M = length(z)

    # Quadrature weighting
    ω = ω_weights(z, θ, σ-1)  # TODO: compute integral for eq 30

    # Discretize the operator. 
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # L_1_minus ≡ L_1 is the only one we use. 

    # Compute the stationary solution.
    stationary_sol = stationary_numerical(params, z) # solve a stationary problem, using numerical solutions
    v_T = stationary_sol.v_tilde
    g_T = stationary_sol.g
    z_hat_T = stationary_sol.z_hat
    Ω_T = stationary_sol.Ω

    # Bundle as before. 
    Ω = get_Ω(Ω_T, δ, T)
    map_z_hat_t = bridge(ℝ, Segment(1, 100))
    map_z_hat_t_inverse = val -> inverse(map_z_hat_t, val) 
    p = @NT(L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ, 
        ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, Ω = Ω, 
        map_z_hat_t = map_z_hat_t, saved_values = SavedValues(Float64, Tuple{Float64,Float64})) #Named tuple for parameters.

    # Dynamic calculations, defined for each time ∈ t.  
    function f!(resid,du,u,p,t)
        @unpack L_1, L_2, z, M, T, μ, υ, Ω, map_z_hat_t = p 

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
    end

    u = [v_T; g_T; map_z_hat_t_inverse(z_hat_T)]
    du = zeros(M+2)
    resid_M2 = zeros(M+2)

    return DAEProblem(f!, resid_M2, u, (T, 0.0), differential_vars = [trues(v_T); false; false], p)

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
    L_tilde_t_derivative = (γ-1)*g_t # default when t = T
    # if (t < T)
    #     throw("L_tilde_t_derivative not defined for t < T") # TODO: remove this and implement derivative fully
    #     forward_index = findlast(x -> x[1] > t, values_future)
    #     t_forward = values_future[forward_index][1]
    #     L_tilde_t_forward = values_future[forward_index][2]
    #     L_tilde_t_derivative = (L_tilde_t_forward - L_tilde_t) / (t_forward - t)
    # end

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
