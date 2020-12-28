# static calculations used by the full model (static and dynamic)
function L_tilde(g, z_hat, Ω, E, S, parameters)
    @unpack N, θ, κ, χ, ζ = parameters
    return Ω * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E / χ)) # (33) (nests the stationary case with E = δ)
end

function L_tilde_x(z_hat, Ω, parameters)
    @unpack N, θ, κ = parameters
    return Ω * (N - 1) * z_hat^(-θ) * κ # (34)
end

L_tilde_E(Ω, E, parameters) = Ω * parameters.ζ * E / parameters.χ # (35) (nests the stationary case with E = δ)
L_tilde_a(Ω, S, parameters) = Ω * parameters.ζ * S # (36)

function S(g, parameters)
    @unpack θ, μ, υ = parameters
    return θ * (g - μ - θ * υ^2 /2) # (32)
end

w(z_bar, parameters) = (parameters.σ - 1)/(parameters.σ) * z_bar # (C.13)

function x(w, parameters)
    @unpack Theta, ζ, η = parameters
    return ζ * (1 - η + η * Theta / w) # (C.14)
end

function λ_ii(z_hat, parameters)
    @unpack N, σ, θ, d = parameters
    return  1/(1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (C.79)
end

c(L_tilde, Ω, z_bar) = (1 - L_tilde)*z_bar # (C.15)
function entry_residual(v_1, Ξ₁, parameters)
    @assert parameters.η ≈ 0.
    Ξ₁*v_1 - parameters.ζ*(1- parameters.χ)/parameters.χ # (56)
end
π_min(L_tilde, z_bar, parameters) = (1 - L_tilde) / ((parameters.σ-1)*z_bar^(parameters.σ-1)) # (38)

function π_rat(z_hat, parameters)
    @unpack θ, σ, d, κ, N = parameters
    return θ/(1 + θ - σ) + (N - 1) * d^(1-σ) * (σ - 1)*z_hat^(σ - 1 - θ)/(1 + θ - σ) # (C.17)
end

function z_bar(z_hat, Ω, parameters)
    @unpack θ, σ, N, d = parameters
    return (Ω*θ/(1 + θ - σ) * (1 + (N - 1)*d^(1 - σ)*z_hat^(σ - 1 - θ)))^(1/(σ-1)) # (C.11)
end

# subset of the above used in the dynamic model
function static_equilibrium(Ξ₁, v_1, g, z_hat, E, Ω, z, parameters)
    @unpack θ, σ, N, d, χ, ζ, κ = parameters
    S_t = S(g, parameters)
    L_tilde_t = L_tilde(g, z_hat, Ω, E, S_t, parameters)
    z_bar_t = z_bar(z_hat, Ω, parameters)
    w_t = w(z_bar_t, parameters)
    π_min_t = (1 - L_tilde_t) / ((σ-1)*z_bar_t^(σ-1))  # (C.20)

    i_vectorized = z .>= log(z_hat) # Vectorized indicator function
    π_t = π_min_t * (1.0.+(N-1)*d^(1-σ)*i_vectorized) - (N-1)*κ*exp.(-(σ-1).*z).*i_vectorized  # (39)
    entry_residual_t = entry_residual(v_1, Ξ₁, parameters)

    return (S = S_t, L_tilde = L_tilde_t, z_bar = z_bar_t, π_min = π_min_t, π = π_t, entry_residual = entry_residual_t,
            w = w_t)
end

consumption_equivalent(U, U_old, parameters) = isapprox(parameters.γ, 1.0) ?  exp(parameters.ρ*(U - U_old)) - 1 : (U/U_old)^(1.0/(1.0 - parameters.γ)) -1# H.23 main paper
ACR(new_λ, old_λ, params) = (-1/params.θ)*log(new_λ/old_λ)
