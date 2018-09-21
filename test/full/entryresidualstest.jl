@testset "entry_residuals Tests" begin 
    # User settings. 
    # Space grid. 
    z_min = 0.0 
    z_max = 5.0 
    M = 1000
    # Experiment settings. 
    d_0 = 5.0
    d_T = 2.3701
    Δ_E = 1e-06
    # Overall parameters. 
    params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053)

    # Construct intermediate objects. 
    z = range(z_min, stop = z_max, length = M)
    params_0 = merge(params, (d = d_0,))
    params_T = merge(params, (d = d_T,))

    # Get (numerical) stationary solution. 
    stationary_0 = stationary_numerical(params_0, z)
    stationary_T = stationary_numerical(params_T, z)

    # Process those 
    Ω_0 = stationary_0.Ω
    Ω_T = stationary_T.Ω
    v_T = stationary_T.v_tilde
    g_T = stationary_T.g 
    z_hat_T = stationary_T.z_hat 

    # Define more interim quantities. 
    T = sqrt(2*(log(Ω_0) - log(Ω_T)) / params.δ) 
    Ω(t) = t < T ? Ω_0 * exp(-params.δ*T*t + params.δ*t*t/2) : Ω_T # Exponential Ω with time smoothing
    E = t -> (log(Ω(t + Δ_E)) - (log(t - Δ_E)))/(2*Δ_E) + params.δ # Log forward differences. 
    
    # Solver settings. 
    tstops = 0:1e-3:T # We don't currently use this anywhere. 
    settings = (z = z, tstops = tstops, Δ_E = Δ_E)

    # Objects for interpolation. 
    Ω_nodes = range(0.0, stop=T, length=30)
    Ω_nodes_interior = Ω_nodes[2:(end-1)]
    entry_residuals_nodes = Ω_nodes

    Ω_interior = map(t -> Ω(t), Ω_nodes_interior)

    # Tests. 
    residuals_interp = entry_residuals(Ω_interior, (Ω_0), stationary_T, T, params_T, settings, Ω_nodes, entry_residuals_nodes).entry_residuals_interpolation
    @test mean(residuals_interp.(entry_residuals_nodes)) ≈ 0.0 atol = 1 # since `Ω` is not in equilibrium, residuals can be high
end 