# Implements the algebraic stationary solution for the full model. Returns the equilibrium quantities (g, Ω, π) determined by equations H.15-H.17.

# Gives us the (full algebraic) stationary solution for a set of params and an initial x.
function stationary_algebraic(params, init_x = defaultiv(params); kwargs...)
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
    c = (1 - L_tilde)*z_bar # (C.74)
    λ_ii = 1/(1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (C.75)
    U_bar = t -> γ == 1 ? (ρ*(log(c) + g*t) + g) / ρ^2 : 1/(1-γ) * (c^(1-γ))/(ρ + (γ-1)*g)  # (C.16)
    (y = c, c = c, λ_ii = λ_ii, U_bar = U_bar)
end

# Gives us the residuals for a point x in state-space and a set of params.
function stationary_algebraic_aux(vals, params)
    # Grab values and intermediate quantities.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals, params)
    g, z_hat, Ω = vals
    # Validate parameters.
    # Calculate and assign residuals.
        big_denom = ν*(θ + ν)*(θ - σ + 1) # (C.18)
        denom_1 = a*(g - r) # (C.18)
        num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # (C.18)
        num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # (C.18)
    return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, # (C.17)
            1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), # (C.18)
            π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))] # (C.19)
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
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # (eq:C.5)
    S = θ * (g - μ - θ * υ^2 /2) # (C.2)
    L_tilde = Ω * ((N-1)*z_hat^(-θ)*κ + (1-η)*ζ*(S + δ/χ)) # (eq:C.7)
    L_tilde_x = Ω * (N - 1) * z_hat^(-θ) * κ # (C.8)
    L_tilde_E = ζ/χ * Ω * δ # (C.9)
    L_tilde_a = ζ * Ω * S # (C.10)
    z_bar = (Ω * (θ/(1 + θ - σ) + (N-1)*(1-F(z_hat))*d^(1-σ)*(z_hat^(-1 + σ)*θ/(1 + θ - σ))))^((σ-1)^(-1)) # (C.11)
    w = σ^(-1)*z_bar # (C,13)
    x = ζ * (1- η + η * Theta / w) # (C.14)
    π_min = (d^(σ-1) * κ)/(z_hat^(σ-1)) # (C.12, inverted)

    return (F = F, r = r, ν = ν, a = a, b = b, S = S, L_tilde = L_tilde, L_tilde_x = L_tilde_x, L_tilde_E = L_tilde_E, L_tilde_a = L_tilde_a,
            z_bar = z_bar, w = w, x = x, π_min = π_min)
end

# Default initial values
defaultiv(params) = [0.02; 18.94; 17.07]

# Numerical method.
function stationary_numerical(params, z, init_x = defaultiv(params); kwargs...)
    # Unpack params and settings.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @assert params.υ > 0 && params.κ > 0 # Parameter validation

    # Discretization objects and quadrature weights.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operators.
    ω = ω_weights(z, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid.

    # Define the system of equations we're solving.
    function stationary_numerical_given_vals(vals)
        g, z_hat, Ω = vals
        @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals([g, z_hat, Ω], params) # Grab static values.
        r_tilde = r - g - 0 # (C.55, and g_w = 0 at steady state)
        ρ_tilde = r_tilde - (σ - 1)*(μ - g + (σ-1)*(υ^2/2)) # (C.40)
        A_T = (ρ_tilde * I - (μ - g + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2) # (46)
        i = z -> z >= log(z_hat) ? 1 : 0 # indicator function for next equation.
        π = z -> π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (33)
        v_tilde = A_T \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (47)

        #=
            System of equations to be solved.
        =#

        # Value-matching condition.
        value_matching = v_tilde[1] - dot(v_tilde, ω) + x # (48) and (C.60)

        # Free-entry condition.
        free_entry = v_tilde[1] - x*(1-χ)/χ # (50) or (C.48) and (C.60)

        # Adoption threshold.
        adoption_threshold = π_min - (1 - L_tilde)/((σ-1)*z_bar^(σ-1)) # (C.19)

        return [value_matching, free_entry, adoption_threshold]
    end

    g_T, z_hat_T, Ω_T = solve_system(stationary_numerical_given_vals, [0.02; 18.94; 17.07])

    # Grab static objects at steady-state.
    staticvalues = staticvals([g_T, z_hat_T, Ω_T], params) # Grab static values.
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvalues
    # Recreate the steady-state objects using the solution in g, z_hat, Ω.
    r_tilde = r - g_T - 0 # (C.55, and g_w = 0 at steady-state)
    ρ_tilde = r_tilde - (σ - 1)*(μ - g_T + (σ-1)*(υ^2/2)) # (C.40)
    A_T = (ρ_tilde * I - (μ-g_T + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2) # (46)
    i = z -> z >= log(z_hat_T) ? 1 : 0 # indicator function for next equation.
    π = z -> π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (33)
    v_tilde = A_T \ π.(z) # discretized system of ODE for v, where v'(T) = 0 (47)

    return merge(staticvalues, merge((g = g_T, z_hat = z_hat_T, Ω = Ω_T, v_tilde = v_tilde), welfare([g_T; z_hat_T; Ω_T], params, staticvalues)))
end
