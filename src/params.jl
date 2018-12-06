parameter_defaults() = @with_kw (ρ = 0.02,
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
                                d = 3.07) # usually set this to 2.5019 for d_T experiment

# some default settings
settings_defaults() = @with_kw (z_grid = unique([range(0., 0.1, length = 400)' range(0.1, 1., length = 400)' range(1., 10, length = 100)']),
                                )
