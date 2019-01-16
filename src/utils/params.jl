parameter_defaults = @with_kw (ρ = 0.02,
                                σ = 3.9896,
                                N = 10,
                                θ = 4.7060,
                                γ = 1.00,
                                κ = 0.0103,
                                ζ = 1.,
                                η = 0.,
                                Theta = 1,
                                χ = 0.4631,
                                υ = 0.0755,
                                μ = 0.,
                                δ = 0.053,
                                d_0 = 2.8,
                                d_T = 2.5019,
                                d = d_T)

# some default settings
settings_defaults = @with_kw (z_max = 5,
                                z = unique([range(0., 0.1, length = 400)' range(0.1, 1., length = 400)' range(1., z_max, length = 100)']),
                                Δ_E = 1e-6,
                                ode_solve_algorithm = CVODE_BDF(),
                                T = 40.0,
                                t = range(0.0, T, length = 10),
                                g = LinearInterpolation(t, stationary_numerical(parameter_defaults(), z).g .+ 0.01*t),
                                E_node_count = 15,
                                entry_residuals_nodes_count = 15,
                                transition_x0 = [-0.91839, -0.797426, -0.61209, -0.500276, -0.319928, -0.235112, -0.235112, -0.121364, -0.100241, -0.0801344, -0.049003, -0.040901, -0.040901, -0.0178076],
                                fifty_node_iv = [-1.00157, -0.848157, -0.821211, -0.821211, -0.821211, -0.748497, -0.633587, -0.527711, -0.498239, -0.498239, -0.498239, -0.498239, -0.3316, -0.3316, -0.3316, -0.3316, -0.3316, -0.281318, -0.281318, -0.281318, -0.281318, -0.281318, -0.241756, -0.230492, -0.168434, -0.168434, -0.168434, -0.168434, -0.105236, -0.103655, -0.103655, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0713765, -0.0713765, -0.0713765, -0.0713765, -0.0343871, -0.0334064, -0.0334064, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373],
                                continuation_x0 = zeros(length(transition_x0)),
                                transition_lb = -ones(length(transition_x0)),
                                transition_ub = zeros(length(transition_x0)),
                                transition_iterations = 1000,
                                transition_penalty_coefficient = 5.0, # coefficient to be used for a penalty function for constraints on increasing E
                                tstops = nothing)

settings_simple_defaults = @with_kw (z = range(0.0, stop = 5.0, length = 100),
                                    T = 100.0,
                                    ode_solve_algorithm = CVODE_BDF(),
                                    iterations = 1000,
                                    t_grid = range(0.0, stop = 100.0, length = length(z)))

parameter_simple_stationary_defaults = @with_kw (μ = 0.0,
    υ = 0.1,
    θ = 2.1,
    r = 0.05,
    ζ = 14.5,
    ξ = 1.0,
    π_tilde = (z->1))

parameter_simple_transition_defaults = @with_kw (μ = 0.0,
    υ = 0.1,
    θ = 2.1,
    r = t -> (0.05 - 1e-02 * (1 - t / 100.0)),
    ζ = 14.5,
    ξ = 1.0,
    π_tilde = (t, z) -> (1 + 1e-02 * (1 - t / 100.0)),
    x = t -> ζ)
