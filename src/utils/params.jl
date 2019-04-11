# Full model objects 
parameter_defaults = @with_kw (ρ = 0.0215,
                                σ = 3.1725, 
                                N = 10,
                                θ = 5.0018,
                                γ = 1.00,
                                κ = 0.0732,
                                ζ = 1.0,
                                η = 0.,
                                Theta = 1,
                                χ = 1/5.9577, 
                                υ = 0.0484 , 
                                μ = -0.0189,
                                δ = 0.02,
                                d = NaN, # must be explicitly set to use this 
                                d_0 = 3.0426,
                                d_T = 1 + (d_0-1)*0.9)

settings_defaults = @with_kw (z_max = 5,
                                z_ex = unique([range(0., 0.1, length = 150)' range(0.1, 1., length = 100)' range(1., 5, length = 50)']),
                                z = z_ex[2:end-1],
                                Δ_E = 1e-6,
                                ode_solve_algorithm = CVODE_BDF(),
                                T = 75.0,
                                t = range(0.0, T, length = 10),
                                weights = [10.0; ones(7)],
                                transition_x0 = [-0.8614905788913748, -0.44926517491157925, -0.283795655883632, -0.0651727069993808, -0.050745349310809824, -0.01213734877329569, -0.005088513022818065, -0.004779208323042127],
                                fifty_node_iv = [-1.00157, -0.848157, -0.821211, -0.821211, -0.821211, -0.748497, -0.633587, -0.527711, -0.498239, -0.498239, -0.498239, -0.498239, -0.3316, -0.3316, -0.3316, -0.3316, -0.3316, -0.281318, -0.281318, -0.281318, -0.281318, -0.281318, -0.241756, -0.230492, -0.168434, -0.168434, -0.168434, -0.168434, -0.105236, -0.103655, -0.103655, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0713765, -0.0713765, -0.0713765, -0.0713765, -0.0343871, -0.0334064, -0.0334064, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373],
                                transition_lb = transition_x0 .- 1e-4,
                                transition_ub = transition_x0 .+ 1e-4,
                                transition_iterations = 10000,
                                transition_penalty_coefficient = 1.0, # coefficient to be used for a penalty function for constraints on increasing E
                                interp = Interpolations.LinearInterpolation,
                                T_U_bar = 50.0,
                                tstops = nothing,
                                impose_E_monotonicity_constraints = true,
                                run_global = true,
                                stationary_x0 = [0.02; 2; .57])

# Simple model objects 
params_simple  = @with_kw (μ = 0.0048, υ = 0.02, θ = 2.1, r = 0.05, ζ = 14.5, ξ = 1.0, π = x -> 1, x = t -> ζ)

settings_simple = @with_kw (z_ex = unique([range(0., stop = 0.1, length = 20)' range(0.1, stop = 1., length = 20)' range(1., stop = 5., length = 20)']), T = 100.0, iterations = 1000, t_grid = range(0., stop = 100., length = length(z_ex)))
