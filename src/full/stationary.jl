# Implements the algebraic stationary solution for the full model. Returns the equilibrium quantities (g, Ω, π) determined by equations H.15-H.17.

# Gives us the (full algebraic) stationary solution for a set of params and an initial x.
function stationary_algebraic(params, init_x = [0.02; 18.94; 17.07]; kwargs...)
    @assert params.υ > 0 && params.κ > 0 # validate params
    g, z_hat, Ω = solve_system(x -> stationary_algebraic_aux(x, params), init_x)
    @assert z_hat > 1 && Ω > 0 && g > 0 # validate solution
    staticvalues = staticvals([g, z_hat, Ω], params)
    return merge(staticvalues, merge((g = g, z_hat = z_hat, Ω = Ω), welfare([g, z_hat, Ω], params, staticvalues))) # calculate quantities and return
end

# Welfare function
function welfare(vals, params, staticvalues)
    g, z_hat, Ω = vals
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals, params)
    c = (1 - L_tilde)*z_bar # (C.15)
    λ_ii = 1/(1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (C.79)
    U_bar = γ == 1 ? (ρ*log(c) + g) / ρ^2 : 1/(1-γ) * (c^(1-γ))/(ρ + (γ-1)*g)  # (C.16)
    return (y = c, c = c, λ_ii = λ_ii, U_bar = U_bar)
end

# Gives us the residuals for a point x in state-space and a set of params.
function stationary_algebraic_aux(vals, params)
    # Grab values and intermediate quantities.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals, params)
    g, z_hat, Ω = vals

    # Calculate and assign residuals.
    big_denom = ν*(θ + ν)*(θ - σ + 1) # (C.19)
    denom_1 = a*(g - r) # (C.19)
    num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # (C.19)
    num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # (C.19)
    return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, # (C.18)
            1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), # (C.19)
            π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))] # (C.20)
end

function staticvals(vals, params)
    # Unpack x.
    g, z_hat, Ω = vals
    # Unpack params.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    # Compute interim quantities.
    F(z) = 1 - z^(-θ) # (C.1)
    r = ρ + γ*g + δ # (C.6)
    ν = (μ-g)/υ^2 + sqrt(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)) # (C.3)
    a = (r - g - (σ - 1)*(μ - g + (σ - 1)*υ^2/2))^(-1) # (C.4)
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # (C.5)
    S = θ * (g - μ - θ * υ^2 /2) # (C.2)
    L_tilde = Ω * ((N-1)*z_hat^(-θ)*κ + (1-η)*ζ*(S + δ/χ)) # (C.7)
    L_tilde_x = Ω * (N - 1) * z_hat^(-θ) * κ # (C.8)
    L_tilde_E = ζ/χ * Ω * δ # (C.9)
    L_tilde_a = ζ * Ω * S # (C.10)
    z_bar = (Ω * (θ/(1 + θ - σ) + (N-1)*(1-F(z_hat))*d^(1-σ)*(z_hat^(-1 + σ)*θ/(1 + θ - σ))))^((σ-1)^(-1)) # (C.11)
    w = (σ-1)/σ*z_bar # (C.13)
    x = ζ * (1- η + η * Theta / w) # (C.14)
    π_min = (d^(σ-1) * κ)/(z_hat^(σ-1)) # (C.12, inverted to express π_min as a function of parameters and z_hat)
    π_rat = (θ + (N-1)*(σ-1)*d^(-θ)*((κ/ζ) * χ/(ρ*(1-χ)))^(1 - θ/(σ - 1)))/(1 + θ - σ) # (C.17)

    return (F = F, r = r, ν = ν, a = a, b = b, S = S, L_tilde = L_tilde, L_tilde_x = L_tilde_x, L_tilde_E = L_tilde_E, L_tilde_a = L_tilde_a,
            z_bar = z_bar, w = w, x = x, π_min = π_min, π_rat = π_rat)
end

# Numerical method checking the stationary solution to the ODE.  In general, use the analtyic when only looking at steady states but this provides guidance on the convergence and accuracy of the finite-difference methods.
function stationary_numerical(params, z_ex, init_x = [0.02; 18.94; 17.07]; kwargs...)
    # Unpack params and settings.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @assert params.υ > 0 && params.κ > 0 # Parameter validation
    z = z_ex[2:end-1] # form a uniform extended grid

    # Discretization objects and quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1_minus = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    ω = ω_weights(z_ex, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid.
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1])) # (24), with ξ = (σ-1)

    # Define the system of equations we're solving.
    function stationary_numerical_given_vals(vals)
        g, z_hat, Ω = vals
        @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals([g, z_hat, Ω], params) # Grab static values.
        r_tilde = r - g - 0 # (C.59, and g_w = 0 at steady state)
        ρ_tilde = r_tilde - (σ - 1)*(μ - g + (σ-1)*(υ^2/2)) # (C.41)
        A_T = ρ_tilde * I - (μ - g + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2 # (52)
        i(z) = z >= log(z_hat) ? 1 : 0 # indicator function for next equation.
        π(z) = π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (39)
        v_tilde = A_T \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (53)

        # System of equations, given the numerically solved ODE
        return [Ξ₁*v_tilde[1] - dot(v_tilde, ω) + ζ, # (54)
                Ξ₁*v_tilde[1] - ζ*(1-χ)/χ, # (56)
                π_min - (1 - L_tilde)/((σ-1)*z_bar^(σ-1))] # (38), from the static equilibrium
    end

    g_T, z_hat_T, Ω_T = solve_system(stationary_numerical_given_vals, [0.02; 18.94; 17.07])

    # Grab static objects at steady-state and recreate the steady-state objects using the g, z_hat, Ω.
    staticvalues = staticvals([g_T, z_hat_T, Ω_T], params)
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvalues
    r_tilde = r - g_T - 0 # (C.59, and g_w = 0 at steady-state)
    ρ_tilde = r_tilde - (σ - 1)*(μ - g_T + (σ-1)*(υ^2/2)) # (C.41)
    A_T = (ρ_tilde * I - (μ-g_T + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2) # (52)
    i(z) = z >= log(z_hat_T) ? 1 : 0 # indicator function for next equation.
    π(z) = π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (39)
    v_tilde = A_T \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (53)

    return merge(staticvalues, merge((g = g_T, z_hat = z_hat_T, Ω = Ω_T, v_tilde = v_tilde), welfare([g_T; z_hat_T; Ω_T], params, staticvalues)))
end
