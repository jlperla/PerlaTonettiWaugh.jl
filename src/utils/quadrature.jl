function ω_weights(z_ex, θ, ξ)
    z_0 = z_ex[1]
    z_bar = z_ex[end]
    Δ_ex = diff(z_ex)
    Δ₋ = [0; Δ_ex] # (19)
    Δ₊ = [Δ_ex; 0] # (20)
    z = z_ex[2:end-1] # i.e., interior of z_ex
    P = length(z)

    ω_bar = 1/2 * (Δ₋ + Δ₊) # (21)
    Ξ₁ = 1/(1 - ξ*Δ₊[1]) # (24)
    Ξₚ = 1/(1 + ξ*Δ₋[end]) # (25)

    ω = similar(z)
    denom = 1 - exp(-θ*z_bar) # (24) denominator (i.e., F(z_bar))

    ω = ω_bar[2:end-1]*θ.*exp.((ξ-θ)*z)/denom # (24)
    ω[1] += ω_bar[1]*Ξ₁*θ*exp((ξ-θ)*z_0)/denom # (24), left endpoint
    ω[end] += ω_bar[end]*Ξₚ*θ*exp((ξ-θ)*z_bar)/denom # (24), right endpoint
    return ω
end
