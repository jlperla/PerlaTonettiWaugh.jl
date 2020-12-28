# Full model objects
parameter_defaults_tests = @with_kw (ρ = 0.0215,
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
                                d_0 = 3.0426,
                                d_T = 2.83834)

settings_defaults_tests = @with_kw (z_ex = unique([range(0., 0.1, length = 150); range(0.1, 1., length = 100); range(1., 5, length = 50)]),
                                T = 75.0,
                                tstops = unique([range(0.0, 10.0, length=41); range(10.0, 20.0, length=21); range(20.0, T, length=56)]),
                                interp = Interpolations.LinearInterpolation,
                                stationary_x0 = default_stationary_x0,
                                pre_shock_times = [-1, -5, -10, -15, -20],
                                fixedpoint_ftol = 1e-8,
                                fixedpoint_show_trace = false,
                                fixedpoint_m = 5,
                                fixedpoint_beta = 1.0,
                                fixedpoint_iterations = 500,
                                fixedpoint_x0 = default_fixedpoint_x0)

parameter_defaults = @with_kw (θ = 4.988976587938262,
                                κ = 0.104196324793307,
                                χ = 1/7.883537319864373,
                                μ = -0.031064624217571,
                                υ = 0.048301140601665,
                                σ = 3.166924135838110,
                                ζ = 1,
                                δ = 0.020,
                                ρ = 0.020338044668517,
                                N = 10,
                                γ = 1.00,
                                η = 0.,
                                Theta = 1,
                                d_0 = 3.022492825462601,
                                d_T = 2.820243542916341)

settings_defaults = @with_kw (z_ex = unique([range(0., 0.1, length = 90); range(0.1, 1., length = 120); range(1., 5, length = 60)]),
                                T = 75.0,
                                tstops = unique([range(0.0, 10.0, length=41); range(10.0, 20.0, length=21); range(20.0, T, length=56)]),
                                interp = Interpolations.LinearInterpolation,
                                stationary_x0 = default_stationary_x0,
                                pre_shock_times = [-1, -5, -10, -15, -20],
                                fixedpoint_ftol = 1e-8,
                                fixedpoint_show_trace = false,
                                fixedpoint_m = 5,
                                fixedpoint_beta = 1.0,
                                fixedpoint_iterations = 500,
                                fixedpoint_x0 = default_fixedpoint_x0)

function default_fixedpoint_x0(parameters, settings)
    if length(settings.tstops) == 116 # HARDCODE THE KNOWN NUMBER OF TSTOPS FOR THE KNOWN SOLUTION
        return [-1.0, -0.988616, -0.977299, -0.96605, -0.954868, -0.943753, -0.932706, -0.921726, -0.910813, -0.899968, -0.88919, -0.878479, -0.867836, -0.85726, -0.846752, -0.836311, -0.825937, -0.81563, -0.805391, -0.79522, -0.785115, -0.775079, -0.765109, -0.755207, -0.745372, -0.735605, -0.725905, -0.716272, -0.706707, -0.69721, -0.687779, -0.678416, -0.669121, -0.659893, -0.650732, -0.641639, -0.632613, -0.623655, -0.614764, -0.605941, -0.588427, -0.571182, -0.554207, -0.5375, -0.521062, -0.504894, -0.488996, -0.473363, -0.458001, -0.442913, -0.42809, -0.413535, -0.399253, -0.385243, -0.371497, -0.358024, -0.344837, -0.331924, -0.319282, -0.306909, -0.2827, -0.2596, -0.237648, -0.216737, -0.196736, -0.177814, -0.160075, -0.143496, -0.12804, -0.113637, -0.100236, -0.0877941, -0.0763549, -0.0659462, -0.0567017, -0.0484165, -0.0408617, -0.0342745, -0.0282847, -0.0234237, -0.0188387, -0.0150505, -0.0118514, -0.00928253, -0.00711861, -0.00539311, -0.00401417, -0.00291141, -0.00208732, -0.00149604, -0.00101103, -0.000673539, -0.000453234, -0.000278941, -0.000168804, -0.0001053, -5.716e-5, -3.03866e-5, -1.73678e-5, -8.03325e-6, -3.65108e-6, -1.99316e-6, -8.91957e-7, -4.8973e-7, -3.65493e-7, -2.81903e-7, -2.38864e-7, -2.06177e-7, -1.73548e-7, -1.40636e-7, -1.07106e-7, -7.23857e-8, -3.64788e-8, 0.0]
    elseif length(settings.tstops) == 88
        return [-1.0, -0.96668, -0.934115, -0.902297, -0.871215, -0.840863, -0.811232, -0.782311, -0.754093, -0.726569, -0.69973, -0.673566, -0.648068, -0.623228, -0.599038, -0.575489, -0.552577, -0.530443, -0.508997, -0.488145, -0.467877, -0.448185, -0.42906, -0.410492, -0.392473, -0.374994, -0.358046, -0.341621, -0.325708, -0.310301, -0.295389, -0.280964, -0.267019, -0.253544, -0.240532, -0.227973, -0.21586, -0.204186, -0.180241, -0.15818, -0.137926, -0.119404, -0.102549, -0.0872978, -0.0735807, -0.0613147, -0.0503832, -0.0407036, -0.0323833, -0.0253242, -0.0195078, -0.0147896, -0.010991, -0.00799515, -0.00567981, -0.00399821, -0.00273242, -0.00183129, -0.00116784, -0.000741846, -0.000442362, -0.000279069, -0.000153617, -8.31759e-5, -4.87385e-5, -2.34855e-5, -1.14183e-5, -6.49374e-6, -3.02169e-6, -1.62968e-6, -1.15147e-6, -8.22077e-7, -6.73772e-7, -5.97719e-7, -5.32128e-7, -4.79283e-7, -4.33382e-7, -3.88373e-7, -3.45597e-7, -3.04862e-7, -2.6164e-7, -2.16769e-7, -1.69966e-7, -1.18744e-7, -6.07212e-8, 0.0]
    else
        return range(-1, 0, length = length(settings.tstops))[2:end-1] #guess linear
    end
end

default_stationary_x0(parameters, settings) = [parameters.ρ * 0.75; 2.0; 0.57]

model_cachename(parameters, settings) = join(hash((
        parameters = parameters,
        settings = merge(settings, (interp = typeof(settings.interp),
                                    stationary_x0 = typeof(settings.stationary_x0),
                                    fixedpoint_x0 = typeof(settings.fixedpoint_x0)))
)))

function load_parameters(filepath::String)
    tmp = load(filepath) |> collect |> first
    return (ρ = tmp.rho, σ = tmp.sigma, θ = tmp.theta, κ = tmp.kappa, χ = tmp.chi, μ = tmp.mu, υ = tmp.upsilon,
            ζ = tmp.zeta, δ = tmp.delta, N = tmp.N, γ = tmp.gamma, η = tmp.eta, Theta = tmp.Theta,
            d_0 = tmp.d_0, d_T = tmp.d_T)
end
