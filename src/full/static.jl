# static calculations used by the full model (static and dynamic)
function L_tilde(g, z_hat, Ω, E, S, params)
    @unpack N, θ, κ, χ, ζ = params
    return Ω * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E / χ)) # (33) (nests the stationary case with E = δ)
end

function L_tilde_x(z_hat, Ω, params)
    @unpack N, θ, κ = params
    return Ω * (N - 1) * z_hat^(-θ) * κ # (34)
end

L_tilde_E(Ω, E, params) = Ω * params.ζ * E / params.χ # (35) (nests the stationary case with E = δ)
L_tilde_a(Ω, S, params) = Ω * params.ζ * S # (36)

function S(g, params)
    @unpack θ, μ, υ = params 
    return θ * (g - μ - θ * υ^2 /2) # (C.2) or (32)
end 

w(z_bar, params) = (params.σ - 1)/(params.σ) * z_bar # (C.13) 

function x(w, params)
    @unpack ζ, Theta, η = params 
    return ζ * (1 - η + η * Theta / w) # (C.14)
end 

function λ_ii(z_hat, params)
    @unpack N, σ, θ, d = params 
    return  1/(1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ)) # (C.79) and (57)
end

c(L_tilde, Ω, z_bar) = (1 - L_tilde)*z_bar # (C.15) and (58)

function z_bar(z_hat, Ω, params)
    @unpack θ, σ, N, d = params 
    return (Ω*θ/(1 + θ - σ) * (1 + (N - 1)*d^(1 - σ)*z_hat^(σ - 1 - θ)))^(1/(σ-1)) # (C.11) and 37 
end 

# subset of the above used in the dynamic model 
function static_equilibrium(Ξ₁, v_1, g, z_hat, E, Ω, z, params)
    @unpack θ, σ, N, d, χ, ζ, κ = params 
    S_t = S(g, params)
    L_tilde_t = L_tilde(g, z_hat, Ω, E, S_t, params)
    z_bar = (Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ)))^(1/(σ-1))  # (37)
    w_t = w(z_bar, params)
    π_min = (1 - L_tilde_t) / ((σ-1)*z_bar^(σ-1))  # (38)

    i_vectorized = z .>= log(z_hat) # Vectorized indicator function
    π = π_min * (1.0.+(N-1)*d^(1-σ)*i_vectorized) - (N-1)*κ*exp.(-(σ-1).*z).*i_vectorized  # (39)
    entry_residual = Ξ₁*v_1 - ζ * (1-χ) / χ  # (56)
    return (S = S_t, L_tilde = L_tilde_t, z_bar = z_bar, π_min = π_min, π = π, entry_residual = entry_residual,
            w = w_t)
end

consumption_equivalent(U, U_old, params) = exp(params.ρ*(U - U_old)) # (C.89) from computational appendix
