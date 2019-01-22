# Diffusion operators with regular grids
function rescaled_diffusionoperators(x::AbstractRange, ξ)
    Δ = step(x)
    P = length(x)

    dl_1 = zeros(P-1)
    d_1 = -ones(P)
    du_1 = ones(P-1)
    d_1[end] = d_1[end] + du_1[end] * (1-ξ*Δ)
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)/Δ

    dl_m1 = -ones(P-1)
    d_m1 = ones(P)
    du_m1 = zeros(P-1)
    d_m1[1] = d_m1[1] + dl_m1[1] * (1+ξ*Δ)
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)/Δ # (A.14) in appendix

    dl_2 = ones(P-1)
    d_2 = -2 * ones(P)
    d_2[1] = -2 + (1+ξ*Δ)
    d_2[end] = -2 + (1-ξ*Δ)
    du_2 = ones(P-1)
    L_2 = Tridiagonal(dl_2, d_2, du_2)/(Δ^2) # (A.15) in appendix

    return (x, L_1_minus, L_1_plus, L_2)
end

# Diffusion operators with irregular grids
function rescaled_diffusionoperators(x::AbstractArray, ξ)
    d = diff(x) # using the first difference as diff from ghost node
    P = length(x)
    Δ_m = zeros(P)
    Δ_m[1] = d[1] 
    Δ_m[2:end] = d
    Δ_p = zeros(P)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    dl_p1 = zeros(P-1)./Δ_m[2:end] 
    d_p1 = -1 * ones(P)./Δ_p
    du_p1 = ones(P-1)./Δ_p[1:end-1]
    d_p1[end] = d_p1[end] + du_p1[end] * (1-ξ*Δ_m[end])
    L_1_plus = Tridiagonal(dl_p1, d_p1, du_p1)

    dl_m1 = -ones(P-1)./Δ_m[2:end]
    d_m1 = ones(P)./Δ_m
    d_m1[1] = d_m1[1] + dl_m1[1] * (1+ξ*Δ_p[1])
    du_m1 = zeros(P-1)./Δ_p[1:end-1]
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)

    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(P-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(P)
    d_2[1] = -2 + (1+ξ * Δ_p[1])
    d_2[end] = -2 + (1-ξ * Δ_m[end])
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(P-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L_2 = Tridiagonal(dl_2, d_2, du_2)
    return (x, L_1_minus, L_1_plus, L_2)
end
