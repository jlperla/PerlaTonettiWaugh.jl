parameter_defaults = @with_kw (ρ = 0.02,
                                σ = 3.05,
                                N = 10,
                                θ = 4.9411,
                                γ = 1.00,
                                κ = 0.1317,
                                ζ = 1.0,
                                η = 0.,
                                Theta = 1,
                                χ = 1/5.2965,
                                υ = 0.0553,
                                μ = -0.0115,
                                δ = 0.05,
                                d = 2.9753)

parameters_old_paper = parameter_defaults(d = 4.0,
                                           θ = 3.1878,
                                           κ = 0.006,
                                           χ = 1.00/2.80,
                                           υ = 0.001,
                                           σ = 3.0,
                                           ζ = 1.00,
                                           δ = 0.001)

# some default settings
settings_defaults = @with_kw (z_max = 5,
                                z_ex = unique([range(0., 0.1, length = 400)' range(0.1, 1., length = 400)' range(1., z_max, length = 100+2)']),
                                z = z_ex[2:end-1],
                                Δ_E = 1e-6,
                                ode_solve_algorithm = CVODE_BDF(),
                                T = 75.0,
                                t = range(0.0, T, length = 10),
                                g = LinearInterpolation(t, stationary_numerical(parameter_defaults(), z).g .+ 0.01*t),
                                weights = [10.0; fill(1.0, 13)],
                                transition_x0 = [-0.8738884967642274, -0.5125465450602664, -0.31163166531947395, -0.1789036699964298, -0.09444764803478085, -0.06676386406358947, -0.03212999344363521, -0.031163943690287882, -0.02126451604463498, -0.006571948963575373, -0.004883226411956144, -0.004613812088943909, -0.003969056948433287, -0.00037064145316761054, -0.0002608582052007029],
                                fifty_node_iv = [-1.00157, -0.848157, -0.821211, -0.821211, -0.821211, -0.748497, -0.633587, -0.527711, -0.498239, -0.498239, -0.498239, -0.498239, -0.3316, -0.3316, -0.3316, -0.3316, -0.3316, -0.281318, -0.281318, -0.281318, -0.281318, -0.281318, -0.241756, -0.230492, -0.168434, -0.168434, -0.168434, -0.168434, -0.105236, -0.103655, -0.103655, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0713765, -0.0713765, -0.0713765, -0.0713765, -0.0343871, -0.0334064, -0.0334064, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373],
                                transition_lb = -ones(length(transition_x0)),
                                transition_ub = zeros(length(transition_x0)),
                                transition_iterations = 1000,
                                transition_penalty_coefficient = 0.0, # coefficient to be used for a penalty function for constraints on increasing E
                                T_U_bar = 50.0,
                                tstops = nothing)

settings_simple_defaults = @with_kw (z_ex = range(0.0, stop = 5.0, length = (100+2)),
                                    z = z_ex[2:end-1],
                                    T = 100.0,
                                    ode_solve_algorithm = CVODE_BDF(),
                                    iterations = 1000,
                                    t_grid = range(0.0, stop = 100.0, length = length(z)))

settings_old_paper_defaults = @with_kw (z_max = 5,
                                z_ex = unique([range(0., 0.1, length = 400)' range(0.1, 1., length = 400)' range(1., 10.0, length = 200+2)']),
                                z = z_ex[2:end-1],
                                Δ_E = 1e-6,
                                ode_solve_algorithm = CVODE_BDF(),
                                T = 120.0,
                                t = range(0.0, T, length = 10),
                                g = LinearInterpolation(t, stationary_numerical(parameter_defaults(), z).g .+ 0.01*t),
                                weights = [10.0; fill(1.0, 13)],
                                transition_x0 = [-0.9802869871313153, -0.7679611162133799, -0.6483822140201239, -0.5709998420691726, -0.4410497194161549, -0.35188633823047205, -0.28134090933192113, -0.22721306548238096, -0.2132657066634307, -0.1802989139615504, -0.1407983331567128, -0.10561616300315106, -0.08546763464126883, -0.058948603687082865, -0.02960294672034148, -0.020649289609280547, -0.013922070758445242, -0.008451708149201357, -0.0039236251615955165],
                                fifty_node_iv = [-1.00157, -0.848157, -0.821211, -0.821211, -0.821211, -0.748497, -0.633587, -0.527711, -0.498239, -0.498239, -0.498239, -0.498239, -0.3316, -0.3316, -0.3316, -0.3316, -0.3316, -0.281318, -0.281318, -0.281318, -0.281318, -0.281318, -0.241756, -0.230492, -0.168434, -0.168434, -0.168434, -0.168434, -0.105236, -0.103655, -0.103655, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0787871, -0.0713765, -0.0713765, -0.0713765, -0.0713765, -0.0343871, -0.0334064, -0.0334064, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373, -0.029373],
                                continuation_x0 = zeros(length(transition_x0)),
                                transition_lb = -ones(length(transition_x0)),
                                transition_ub = zeros(length(transition_x0)),
                                transition_iterations = 1000,
                                transition_penalty_coefficient = 0.0, # coefficient to be used for a penalty function for constraints on increasing E
                                T_U_bar = 50.0,
                                tstops = nothing)

parameter_simple_stationary_defaults = @with_kw (μ = 0.0,
    υ = 0.1,
    θ = 2.1,
    r = 0.05,
    ζ = 14.5,
    ξ = 1.0)

parameter_simple_transition_defaults = @with_kw (μ = 0.0,
    υ = 0.1,
    θ = 2.1,
    r = t -> (0.05 - 1e-02 * (1 - t / 100.0)),
    ζ = 14.5,
    ξ = 1.0,
    x = t -> ζ)
