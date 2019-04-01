# New method.
function ω_weights(z_ex, θ, ξ)
    # preliminaries
    # grid objects
        z_0 = z_ex[1]
        z_bar = z_ex[end]
        z = z_ex[2:end-1] # i.e., interior of z_ex
        d = diff(z_ex)
        Δ₋ = [0; d] # (A.3)
        Δ₊ = [d; 0] # (A.4)
        P = length(z)
    # ω objects
        ω = zeros(eltype(z), P) # object to return
        ω_bar_vec = 1/2 * (Δ₋ + Δ₊) # (B.19)
        ω_bar(i) = ω_bar_vec[i+1] # so we can match the notation of (B.21) and write ω_bar(0), ω_bar(1), ..., ω_bar(P+1)

    # fill ω
    for i = 1:P
        if i == 1
            Ξ₁ = 1/(1 - ξ*(z[1] - z_0)) # (A.11)
            ω[1] = ω_bar(0) * Ξ₁ * (θ*exp( (ξ-θ)*z_0 ))/(1 - exp(-θ * z_bar)) + ω_bar(1)*(exp( (ξ-θ)*z[1] ))/(1 - exp(-θ * z_bar)) # (24), left endpoint
        elseif i == P
            Ξₚ = 1/(1 + ξ*(z_ex[end] - z_ex[end-1])) # (A.12)
            ω[P] = ω_bar(P) * (θ*exp( (ξ-θ)*z[P] ))/(1 - exp(-θ * z_bar)) + ω_bar(P+1)*Ξₚ*(θ*exp( (ξ-θ)*z_bar ))/(1 - exp(-θ * z_bar))
        else
            ω[i] = ω_bar(i) * (θ * exp( (ξ - θ)*z[i] ))/(1 - exp(-θ * z_bar))
        end
    end
    return ω
end
