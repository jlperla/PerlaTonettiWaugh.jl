
#These would be generated automatically by DiffEqOperators.jl
#These are the differential operators for positive drift 0upwind finite differences
function diffusionoperators(x_min, x_max, M)
    x = linspace(x_min, x_max, M)
    Δ = x[2]-x[1]

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
    return (x, BandedMatrix(L_1_minus, (1,1)), BandedMatrix(L_1_plus, (1, 1)), BandedMatrix(L_2, (1, 1))) #The (1,1) are the off-diagonal bandwidths
end
