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
                                d_0 = 3.0426,
                                d_T = 1 + (d_0-1)*0.9)

settings_defaults = @with_kw (z_ex = unique([range(0., 0.1, length = 150)' range(0.1, 1., length = 100)' range(1., 5, length = 50)']),
                                T = 75.0,
                                tstops = unique([collect(0.0:0.25:10); collect(10: 0.5: 20.); collect(20.0:T)]),
                                interp = Interpolations.LinearInterpolation,
                                stationary_x0 = default_stationary_x0,
                                pre_shock_times = [-1, -5, -10, -15, -20],
                                fixedpoint_ftol = 1e-8,
                                fixedpoint_show_trace = false,
                                fixedpoint_m = 5,
                                fixedpoint_beta = 1.0,
                                fixedpoint_iterations = 500,
                                fixedpoint_x0 = default_fixedpoint_x0)

default_fixedpoint_x0(parameters, settings) = range(-1, 0, length = length(settings.tstops))[2:end-1]
default_stationary_x0(parameters, settings) = [parameters.ρ * 0.75; 2.0; 0.57]

model_cachename(parameters, settings) = join(hash((
        parameters = parameters, 
        settings = merge(settings, (interp = typeof(settings.interp), 
                                    stationary_x0 = typeof(settings.stationary_x0), 
                                    fixedpoint_x0 = typeof(settings.fixedpoint_x0)))
)))