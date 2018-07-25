# Method for regular grids. 
function diffusionoperators(z::Range)
    Δ = step(z)
    M = length(z)

    dl_1 = zeros(M-1)
    d_1 = -1 * ones(M)
    d_1[end] = 0
    du_1 = ones(M-1)
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)/Δ

    dl_m1 = -ones(M-1)/Δ
    d_m1 = ones(M)/Δ
    d_m1[1] = 0 
    du_m1 = zeros(M-1)/Δ
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)

    dl_2 = ones(M-1)
    d_2 = -2 * ones(M)
    d_2[1] = -1
    d_2[end] = -1
    du_2 = ones(M-1)
    L_2 = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)

    #BandedMatrix are much faster, probably because of better specializations in the composition
    return (z, BandedMatrix(L_1_minus, (1,1)), BandedMatrix(L_1_plus, (1, 1)), BandedMatrix(L_2, (1, 1))) #The (1,1) are the off-diagonal bandwidths
end 

# Method for irregular grids. 
function diffusionoperators(z::AbstractArray)
    d = diff(z)
    M = length(z)
    Δ_m = zeros(M)
    Δ_m[1] = d[1] # using the first difference as diff from ghost node
    Δ_m[2:end] = d
    Δ_p = zeros(M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    # I think l denote Z line, u dentoe X line
    dl_p1 = zeros(M-1)./Δ_m[2:end] # This is flipped delta vector as the tridiagnal matrix is flipped? not sure
    d_p1 = -1 * ones(M)./Δ_p
    d_p1[end] = 0
    du_p1 = ones(M-1)./Δ_p[1:end-1]
    L_1_plus = Tridiagonal(dl_p1, d_p1, du_p1)

    dl_m1 = -ones(M-1)./Δ_m[2:end]
    d_m1 = ones(M)./Δ_m
    d_m1[1] = 0
    du_m1 = zeros(M-1)./Δ_p[1:end-1]
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)

    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(M-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(M)
    d_2[1] = -1
    d_2[end] = -1
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(M-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L_2 = Tridiagonal(dl_2, d_2, du_2)

    #BandedMatrix are much faster, probably because of better specializations in the composition
    return (z, BandedMatrix(L_1_minus, (1, 1)),BandedMatrix(L_1_plus, (1, 1)), BandedMatrix(L_2, (1, 1))) #The (1,1) are the off-diagonal bandwidths
end 

