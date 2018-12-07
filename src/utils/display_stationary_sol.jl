# auxil. function that extracts and displays welfare information from steady state solutions from full models
function display_stationary_sol(stationary_sol)
    @unpack g, z_hat, Ω, π_bar_agg, y, c, λ_ii, U_bar, L_tilde, z_bar, w, x, π_min, r, a, b, S = stationary_sol
    @show g
    @show z_hat
    @show Ω
    @show π_bar_agg
    @show y
    @show c
    @show λ_ii
    @show U_bar(0.0)
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