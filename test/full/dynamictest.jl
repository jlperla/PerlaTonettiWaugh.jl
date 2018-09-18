# User settings. 
  # Time and space grid. 
    z_min = 0.0 
    z_max = 5.0 
    M = 1000
    T = 20.0 
  # Experiment settings. 
    d_0 = 5.0
    d_T = 2.3701
    Δ_E = 1e-06
  # Overall parameters. 
    params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053)
  # Solver settings. 
    tstops = 0:1e-3:T 

# Construct intermediate objects. 
  z = range(z_min, stop = z_max, length = M)
  settings = (z = z, tstops = tstops, Δ_E = Δ_E)
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
  Ω = t -> Ω_T # This is constant. 
  E = t -> (log(Ω(t + \Delta) ))