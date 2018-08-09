# Implements the algebraic stationary solution for the full model. Returns the equilibrium quantities (g, Ω, π) determined by equations H.15-H.17.


# Gives us the (full algebraic) stationary solution for a set of params and an initial x. 
stationary_algebraic_full(init_x, params; kwargs...) = _staticvalues(nlsolve(x -> _fullstationary(x, params), init_x; inplace = false, kwargs...).zero, params)

# Gives us the residuals for a point x in state-space and a set of params. 
function _fullstationary(x, params)    
    # Grab values and intermediate quantities. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack g, z_hat, Ω, F, r, ν, a, b, S, L_tilde, z_bar, w, x, π = _staticvalues(x, params)
    # Validate parameters. 
    @assert υ > 0 
    @assert κ > 0
    @assert z_hat > 1 && Ω > 0 && g > 0
    # Calculate and assign residuals. 
    G = zeros(3)
    G[1] = x/π - a*(χ/(1-χ))*(σ + ν - 1)/ν # H.15
        # Calculate H.16 intermediates. 
        big_denom = ν*(θ + ν)*(θ - σ + 1) # Part of H.16
        denom_1 = a*(g - r) # Part of H.16
        num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # Part of H.16
        num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # Part of H.16
    G[2] = 1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν) # H.16
    G[3] = π - (1- L_tilde)/((σ -1)*z_bar^(σ-1)) # H.17
    return G
end

function _staticvalues(x, params)
    # Unpack x. 
    g = x[1]
    z_hat = x[2]
    Ω = x[3]
    # Unpack params. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    # Compute interim quantities.
    F(z) = 1 - z^(-θ) # H.1 
    r = ρ + γ*g + δ # H.6 
    ν = (μ-g)/υ^2 + sqrt(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)) # H.3
    a = inv(r - g - (σ - 1)*(μ - g + (σ - 1)*υ^2/2)) # H.4 
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # H.5 
    S = θ * (g - μ - θ * υ^2 /2) # H.2   
    L_tilde = Ω * ((N-1)*(1-F(z_hat))*κ + (1-η)*ζ*(S + δ/χ)) # H.7
    z_bar = (Ω * (θ/(1 + θ - σ) + (N-1)*(1-F(z_hat))*d^(1-σ)*(z_hat^(-1 + σ)*θ/(1 + θ - σ))))^(inv(σ-1)) # H.8
    w = inv(σ)*z_bar # H.10
    x = ζ * (1- η + η * Theta / w) # H.11
    π = (d^(σ-1) * κ)/(z_hat^(σ-1)) # Inversion of H.9    
    return @NT(g = g, z_hat = z_hat, Ω = Ω, F = F, r = r, ν = ν, a = a, b = b, S = S, L_tilde = L_tilde, z_bar = z_bar, w = w, x = x, π = π)
end 