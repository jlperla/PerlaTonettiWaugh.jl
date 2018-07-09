#These would be generated automatically by DiffEqOperators.jl
#These are the differential operators for positive drift 0upwind finite differences
function irregulardiffusionoperators(x, M)

    d = diff(x);
    Δ_m = zeros(M);
    Δ_m[1] = d[1]; # using the first difference as diff from ghost node
    Δ_m[2:end] = d;
    Δ_p = zeros(M);
    Δ_p[end] = d[end];
    Δ_p[1:end-1] = d;

    # I think l denote Z line, u dentoe X line
    dl_1 = zeros(M-1)./Δ_p[end:-1:2] # This is flipped delta vector as the tridiagnal matrix is flipped? not sure
    d_1 = -1 * ones(M)./Δ_p[end:-1:1]
    d_1[end] = 0
    du_1 = ones(M-1)./Δ_p[end:-1:2]
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)

    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(M-1)./(Δ_m[end:-1:2].*Δ[end:-1:2])
    d_2 = -2 * ones(M)
    d_2[1] = -1
    d_2[end] = -1
    d_2 = d_2./(Δ_p[end:-1:1].*Δ_m[end:-1:1])
    du_2 = 2*ones(M-1)./(Δ_p[end:-1:2].*Δ[end:-1:2])
    L_2 = Tridiagonal(dl_2, d_2, du_2)

    #BandedMatrix are much faster, probably because of better specializations in the composition
    return (x, BandedMatrix(L_1_plus, (1, 1)), BandedMatrix(L_2, (1, 1))) #The (1,1) are the off-diagonal bandwidths
end