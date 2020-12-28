function total_derivative(params_baseline, ϵ = 0.000001, settings = settings_defaults())
    @unpack ρ, d, γ = params_baseline
    sol_baseline = stationary_algebraic(params_baseline, settings)
    if γ == 1
        U_1 = 1/(ρ*sol_baseline.c)
        U_2 = (1/ρ^2)
    else
        U_1 = sol_baseline.c^(-γ) / (ρ + (γ-1)* sol_baseline.g)
        U_2 = sol_baseline.c^(1-γ) * (ρ + (γ-1)* sol_baseline.g)^(-2)
    end

    d_counterfactual = d + ϵ
    params_counterfactual = merge(params_baseline, (d = d_counterfactual,))
    sol_counterfactual = stationary_algebraic(params_counterfactual, settings)

	ACR = (-1/params_baseline.θ)*log(sol_counterfactual.λ_ii/sol_baseline.λ_ii)

    sol_fc_d = steady_state_from_g(sol_baseline.g, sol_baseline.z_hat, sol_baseline.Ω, params_counterfactual, settings)
    sol_fc_g = steady_state_from_g(sol_baseline.g + ϵ, sol_baseline.z_hat, sol_baseline.Ω, params_baseline, settings)
    sol_fc_zhat = steady_state_from_g(sol_baseline.g, sol_baseline.z_hat + ϵ, sol_baseline.Ω, params_baseline, settings)
    sol_fc_Omega = steady_state_from_g(sol_baseline.g, sol_baseline.z_hat, sol_baseline.Ω + ϵ, params_baseline, settings)

    partial_fc_d = (sol_fc_d.c - sol_baseline.c)/(d_counterfactual - params_baseline.d)
    partial_fc_g = (sol_fc_g.c - sol_baseline.c)/(sol_fc_g.g - sol_baseline.g)
    partial_fc_zhat = (sol_fc_zhat.c - sol_baseline.c)/(sol_fc_zhat.z_hat - sol_baseline.z_hat)
    partial_fc_Omega = (sol_fc_Omega.c - sol_baseline.c)/(sol_fc_Omega.Ω - sol_baseline.Ω)

    partial_Omega_d = (sol_counterfactual.Ω - sol_baseline.Ω)/(d_counterfactual - params_baseline.d)
    partial_zhat_d = (sol_counterfactual.z_hat - sol_baseline.z_hat)/(d_counterfactual - params_baseline.d)
    partial_g_d = (sol_counterfactual.g - sol_baseline.g)/(d_counterfactual - params_baseline.d)

    d_U_d_total = (sol_counterfactual.U_bar - sol_baseline.U_bar)/(d_counterfactual - params_baseline.d)
    total_decomp = U_1*(partial_fc_Omega*partial_Omega_d + partial_fc_zhat*partial_zhat_d + partial_fc_g*partial_g_d + partial_fc_d) + U_2*partial_g_d
    check = total_decomp/d_U_d_total
	total_decomp_CE = ρ*total_decomp

    decomp_fc_d = partial_fc_d
    decomp_fc_Omega_Omega_d = partial_fc_Omega*partial_Omega_d
    decomp_fc_zhat_zhat_d = partial_fc_zhat*partial_zhat_d
    decomp_fc_g_g_d = partial_fc_g*partial_g_d
    decomp_g_d = partial_g_d

    U1_partial_fc_d_frac = U_1*partial_fc_d/total_decomp
    U1_decomp_fc_Omega_Omega_d_frac = U_1*decomp_fc_Omega_Omega_d/total_decomp
	U1_decomp_fc_zhat_zhat_d_frac = U_1*decomp_fc_zhat_zhat_d/total_decomp
	U1_decomp_fc_g_g_d_frac = U_1*decomp_fc_g_g_d/total_decomp
	U2_decomp_g_d_frac = U_2*decomp_g_d/total_decomp

    planner_0_g_CE = ρ*(U_1*partial_fc_g +U_2)
    planner_0_Omega_CE = ρ*U_1*partial_fc_Omega
    planner_0_zhat_CE = ρ*U_1*partial_fc_zhat
	MU_g_CE = (1/ρ)
	semi_elasticity_g_d = params_baseline.d*partial_g_d
	U1_partial_fc_d_ACRfrac = U1_partial_fc_d_ACRfrac = ρ*U_1*partial_fc_d*ϵ/ACR
	d_T = 1 + 0.9*(params_baseline.d - 1)
	d_d = d_T - params_baseline.d
	CE_extrap10p = total_decomp_CE * d_d
	g_extrap10p = partial_g_d * d_d

    return (check = check, d_U_d_total = d_U_d_total, U_1 = U_1, U_2 = U_2,
										ACR = ACR,
										∂_fc_d = partial_fc_d,
										∂_fc_g = partial_fc_g,
										∂_fc_zhat = partial_fc_zhat,
										∂_fc_Omega = partial_fc_Omega,
										∂_g_d = partial_g_d,
										∂_zhat_d = partial_zhat_d,
										∂_Omega_d = partial_Omega_d,
										total_decomp = total_decomp,
										decomp_fc_d = decomp_fc_d,
										decomp_fc_Omega_Omega_d = decomp_fc_Omega_Omega_d ,
										decomp_fc_zhat_zhat_d = decomp_fc_zhat_zhat_d,
										decomp_fc_g_g_d = decomp_fc_g_g_d,
										decomp_g_d = decomp_g_d,
										U1_partial_fc_d_frac = U1_partial_fc_d_frac,
										U1_decomp_fc_Omega_Omega_d_frac = U1_decomp_fc_Omega_Omega_d_frac,
										U1_decomp_fc_zhat_zhat_d_frac = U1_decomp_fc_zhat_zhat_d_frac,
										U1_decomp_fc_g_g_d_frac = U1_decomp_fc_g_g_d_frac,
										U2_decomp_g_d_frac = U2_decomp_g_d_frac,
										planner_0_g_CE = planner_0_g_CE,
										planner_0_Omega_CE = planner_0_Omega_CE,
										planner_0_zhat_CE = planner_0_zhat_CE,
										MU_g_CE = MU_g_CE,
										semi_elasticity_g_d = semi_elasticity_g_d,
										U1_partial_fc_d_ACRfrac = U1_partial_fc_d_ACRfrac,
										d_T = d_T,
										d_d = d_d,
										CE_extrap10p = CE_extrap10p,
										g_extrap10p = g_extrap10p)
end

function steady_state_from_g(g, z_hat, Ω, parameters, settings)
    @assert parameters.υ > 0 && parameters.κ > 0
    @assert z_hat > 1 && Ω > 0 && g > 0 # input validation still required
    staticvalues = staticvals([g, z_hat, Ω], parameters)
    return merge(staticvalues, merge((g = g, z_hat = z_hat, Ω = Ω,), welfare([g, z_hat, Ω], parameters)))
end

function stationary_algebraic(parameters, settings)
    x0 = settings.stationary_x0(parameters, settings)
    @assert parameters.υ > 0 && parameters.κ > 0
    sol = nlsolve(x -> stationary_algebraic_aux(x, parameters), x0, inplace = false) # uses distorted parameters
    converged(sol) || error("Solver didn't converge.")
    g, z_hat, Ω = sol.zero
    @assert z_hat > 1 && Ω > 0 && g > 0
    staticvalues = staticvals([g, z_hat, Ω], parameters) # uses undistorted parameters (except x)

    # calculate value function using differential objects based on the grid
    @unpack z_ex = settings
    z = z_ex[2:end-1]
    @unpack θ, σ = parameters
    bc = (Mixed(ξ = σ-1), Mixed(ξ = σ-1)) # boundary conditions for differential operators
    L_1_minus = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    ω = ω_weights(z_ex, θ, σ-1) # get quadrature weights for the distribution on the rescaled grid.
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1])) # (24), with ξ = (σ-1)
    p = (L_1_minus = L_1_minus, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, z = z)
    v_tilde = stationary_numerical_given_vals([g, z_hat, Ω], p, parameters, settings).v_tilde

    return merge(staticvalues, merge((g = g, z_hat = z_hat, Ω = Ω, v_tilde = v_tilde), welfare([g, z_hat, Ω], parameters))) # calculate quantities and return
end

function welfare(vals, parameters)
    g, z_hat, Ω = vals
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, π_min = staticvals(vals, parameters)
    c_T = c(L_tilde, Ω, z_bar)
    return (y = c_T, # y = c by (C.78)
            c = c_T,
            λ_ii = λ_ii(z_hat, parameters),
            U_bar = γ == 1 ? (ρ*log(c_T) + g) / ρ^2 : 1/(1-γ) * (c_T^(1-γ))/(ρ + (γ-1)*g)) # (C.16)
end

# Kernel for the algebraic stationary.
function stationary_algebraic_aux(vals, parameters)
    @unpack ρ, σ, N, θ, γ, d, κ, η, Theta, υ, μ, δ, χ, ζ = parameters
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals, parameters)
    g, z_hat, Ω = vals

    big_denom = ν*(θ + ν)*(θ - σ + 1) # (C.19)
    denom_1 = a*(g - r) # (C.19)
    num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # (C.19)
    num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # (C.19)
    return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, # (C.18)
            1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), # (C.19)
            π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))] # (C.20)
end

function staticvals(vals, parameters)
    g, z_hat, Ω = vals
    @assert parameters.η ≈ 0  # we think that the L_tilde formulas are wrong for η != 0
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters

    F(z) = 1 - z^(-θ) # (C.1)
    r = ρ + γ*g + δ # (C.6)
    ν = (μ-g)/υ^2 + sqrt(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)) # (C.3)
    a = (r - g - (σ - 1)*(μ - g + (σ - 1)*υ^2/2))^(-1) # (C.4)
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # (C.5)
    S_T = S(g, parameters) # calls like this rely on utils/static.jl
    L_tilde_T = L_tilde(g, z_hat, Ω, δ, S_T, parameters)
    L_tilde_x_T = L_tilde_x(z_hat, Ω, parameters)
    L_tilde_E_T = L_tilde_E(Ω, δ, parameters)
    L_tilde_a_T = L_tilde_a(Ω, S_T, parameters)
    z_bar_T = z_bar(z_hat, Ω, parameters)
    w_T = w(z_bar_T, parameters)
    x_T = x(w_T, parameters)
    π_min_T = (d^(σ-1) * κ)/(z_hat^(σ-1)) # (C.12, inverted to express π_min as a function of parameters and z_hat) This is where the export threshold is used.  pi_min(L_tilde, z_bar)  is already used in the system of equations
    π_rat_T = π_rat(z_hat, parameters)

    return (F = F, r = r, ν = ν, a = a, b = b, S = S_T, L_tilde = L_tilde_T, L_tilde_x = L_tilde_x_T, L_tilde_E = L_tilde_E_T, L_tilde_a = L_tilde_a_T,
            z_bar = z_bar_T, w = w_T, x = x_T, π_min = π_min_T, π_rat = π_rat_T)
end

# Kernel of the stationary numerical.
function stationary_numerical_given_vals(vals, p, parameters, settings)
    g, z_hat, Ω = vals
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, π_min = staticvals([g, z_hat, Ω], parameters)
    @unpack ρ, σ, N, θ, γ, d, κ, η, Theta, ζ, χ, υ, μ, δ = parameters
    @unpack Ξ₁, L_1_minus, L_2, ω, z = p
    r_tilde = r - g - 0 # (C.59, and g_w = 0 at steady state)
    ρ_tilde = r_tilde - (σ - 1)*(μ - g + (σ-1)*(υ^2/2)) # (C.41)
    A_T = ρ_tilde * I - (μ - g + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2 # (52)
    i(z) = z >= log(z_hat) ? 1 : 0 # indicator function for next equation.
    π(z) = π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (39)
    v_tilde = A_T \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (53)

    # System of equations, given the numerically solved ODE
    return (residuals = [Ξ₁*v_tilde[1] - dot(v_tilde, ω) + ζ, # (54)
            Ξ₁*v_tilde[1] - ζ*(1-χ)/χ, # (56)
            π_min - (1 - L_tilde)/((σ-1)*z_bar^(σ-1))], # (38)
            v_tilde = v_tilde)
end

# Solve the ODE by root-finding. Should use the analytic most of the time, but this is a useful check.
function stationary_numerical(parameters, settings)
    @unpack z_ex = settings
    x0 = settings.stationary_x0(parameters, settings)
    @unpack ρ, σ, N, θ, γ, d, κ, η, Theta, υ, μ, δ = parameters
    @assert υ > 0 && κ > 0 # Parameter validation
    z = z_ex[2:end-1] # form a uniform extended grid

    # Differential objects we need for the residuals
    bc = (Mixed(ξ = σ-1), Mixed(ξ = σ-1)) # boundary conditions for differential operators
    L_1_minus = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    ω = ω_weights(z_ex, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid.
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1])) # (24), with ξ = (σ-1)

    p = (Ξ₁ = Ξ₁, L_1_minus = L_1_minus, L_2 = L_2, ω = ω, z = z)
    sol = nlsolve(vals -> stationary_numerical_given_vals(vals, p, parameters, settings).residuals, x0, inplace = false)
    converged(sol) || error("Solver didn't converge.")
    g_T, z_hat_T, Ω_T = sol.zero

    staticvalues = staticvals([g_T, z_hat_T, Ω_T], parameters)
    v_tilde = stationary_numerical_given_vals([g_T, z_hat_T, Ω_T], p, parameters, settings).v_tilde
    return merge(staticvalues, merge((g = g_T, z_hat = z_hat_T, Ω = Ω_T, v_tilde = v_tilde), welfare([g_T; z_hat_T; Ω_T], parameters)))
end

# Utilities for comparison, display, etc.
function compare_steady_states(parameters, settings; verbose = false, algebraic = false)

    if algebraic
        stationary = stationary_algebraic
    else
        stationary = stationary_numerical
    end

    parameters_0 = merge(parameters, (d = parameters.d_0,))
    parameters_T = merge(parameters, (d = parameters.d_T,))
    stationary_0 = stationary(parameters_0, settings)
    stationary_T = stationary(parameters_T, settings)
    change_welfare = 100*consumption_equivalent(stationary_T.U_bar, stationary_0.U_bar, parameters)
    change_trade = (1-stationary_T.λ_ii) - (1-stationary_0.λ_ii)

    if verbose
        println("SS to SS welfare gain: $(change_welfare)")
        println("Change in Trade: $(change_trade)")
        println("Growth Rates across SS: $(100*stationary_T.g) (time T) vs $(100*stationary_0.g) (time 0)")
    end

    return (stationary_0 = stationary_0,
            stationary_T = stationary_T,
            change_welfare = change_welfare,
            change_trade = change_trade)
end
