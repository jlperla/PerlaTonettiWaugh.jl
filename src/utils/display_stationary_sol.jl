# auxil. function that extracts and displays welfare information from steady state solutions from full models
function display_stationary_sol(stationary_sol)
    @unpack g, z_hat, Ω, y, c, λ_ii, L_tilde, U_bar, z_bar, w, x, π_min, r, a, b, S = stationary_sol
    @show g
    @show z_hat
    @show Ω
    @show y
    @show c
    @show U_bar
    @show λ_ii
    @show L_tilde
    @show z_bar
    @show w
    @show x
    @show π_min
    @show r
    @show a
    @show b
    @show S
end;
