# Simple model objects
parameters_simple  = @with_kw (μ = 0.0048, 
                                υ = 0.02, 
                                θ = 2.1, 
                                r = t -> 0.05, 
                                ζ = 14.5,
                                ξ = 1.0, π = (t, x) -> 1, 
                                x = t -> ζ)

settings_simple = @with_kw (z_ex = unique([range(0., 0.1, length = 20)' range(0.1, 1., length = 20)' range(1., 5., length = 20)']),
                            T = 100.0, iterations = 1000, 
                            ts = range(0., 100., length = length(z_ex)), 
                            x0 = [0.1], 
                            ode_solve_algorithm = IDA())
