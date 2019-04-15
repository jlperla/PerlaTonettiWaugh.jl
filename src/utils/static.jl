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
F(params) = z -> 1 - z^(-params.θ) # (C.1)

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
