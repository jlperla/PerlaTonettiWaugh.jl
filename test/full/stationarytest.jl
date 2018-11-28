# User settings. 
    # State grids. 
    grids = [   range(0.0, 5.0, length = 500),
                range(0.0, 7.0, length = 700),
                unique([range(0.0, 5.0, length = 500)' range(5.0, 7.0, length = 200)']), # irregular 
                unique([range(0.0, 5.0, length = 500)' range(0.0, 5.0, length = 330)']),
            ]
    # Overall parameters. 
    params = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1,                    η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) 
    baseline = params()
    # Solver settings. 
    initial_values = [  [0.25, 3.0, 1.0],
                        [0.0190, 1.434969, 1.06517] # ~ equilibrium value for these parameters.
                     ]

# Execute tests. 
@testset "Tests with Vanilla Parameters" begin 
    # Compute algebraic solutions. 
    algebraic_sols = stationary_algebraic.(Ref(baseline), initial_values)
    # Test algebraic equilibrium quantities. 
    # Growth rate tests. 
    algebraic_gs = (x -> x.g).(algebraic_sols)
    @test var(algebraic_gs) < 1e-10 # Tests that the solutions are similar to one another. 
    @test all(algebraic_gs .≈ 0.01900455065415125) # Tests proximity to true value. 
    # Ω tests. 
    algebraic_Ωs = (x -> x.Ω).(algebraic_sols)
    @test var(algebraic_Ωs) < 1e-10
    @test all(algebraic_Ωs .≈ 1.06517755)
    # z_hat tests. 
    algebraic_zs = (x -> x.z_hat).(algebraic_sols)
    @test var(algebraic_zs) < 1e-10
    @test all(algebraic_zs .≈ 1.434969541725385)

    # Compute numerical solutions.
    numerical_sols = stationary_numerical.(Ref(baseline), grids)
    # Test numerical equilibrium quantities. 
    numerical_gs = (x -> x.g).(numerical_sols)
    @test_broken var(numerical_gs) < 1e-10 # Tests that the solutions are similar to one another. 
    @test_broken all(numerical_gs .≈ 0.01900455065415125) # Tests proximity to true value. 
    # Ω tests. 
    numerical_Ωs = (x -> x.Ω).(numerical_sols)
    @test_broken var(numerical_Ωs) < 1e-10
    @test_broken all(numerical_Ωs .≈ 1.06517755)
    # z_hat tests. 
    numerical_zs = (x -> x.z_hat).(numerical_sols)
    @test_broken var(numerical_zs) < 1e-10
    @test_broken all(numerical_zs .≈ 1.434969541725385)

    # Numerical residuals tests. 
    for i in 1:length(numerical_sols)
        # Get values to test. 
        grid = grids[i]
        sol = numerical_sols[i]
        # Compute interim quantities. 
        ω = ω_weights(grid, baseline.θ, baseline.σ-1)
        value_matching = sol.v_tilde[1] - dot(sol.v_tilde, ω) + sol.x
        free_entry = sol.v_tilde[1] - sol.x*(1-baseline.χ)/baseline.χ
        # Test. 
        @test value_matching < 1e-8
        @test_broken free_entry < 1e-8
    end     
end 

@testset "Tests with Perturbed Parameters" begin 
    # Generate new parameters. 
    newparams = params(σ = 4.0)
    # Compute algebraic solutions. 
    algebraic_sols = stationary_algebraic.(Ref(newparams), initial_values)
    # Test consistency. 
    algebraic_gs = (x -> x.g).(algebraic_sols)
    algebraic_Ωs = (x -> x.Ω).(algebraic_sols)
    algebraic_zs = (x -> x.z_hat).(algebraic_sols)
    @test var(algebraic_gs) < 1e-10 # Tests that the solutions are similar to one another. 
    @test var(algebraic_Ωs) < 1e-10
    @test var(algebraic_zs) < 1e-10
    
    # Compute numerical solutions. 
    numerical_sols = stationary_numerical.(Ref(newparams), grids)
    # Test consistency. 
    numerical_gs = (x -> x.g).(numerical_sols)
    numerical_Ωs = (x -> x.Ω).(numerical_sols)
    numerical_zs = (x -> x.z_hat).(numerical_sols)
    @test_broken var(numerical_gs) < 1e-10 # Tests that the solutions are similar to one another. 
    @test_broken var(numerical_Ωs) < 1e-10
    @test_broken var(numerical_zs) < 1e-10

    # Numerical residuals tests. 
    for i in 1:length(numerical_sols)
        # Get values to test. 
        grid = grids[i]
        sol = numerical_sols[i]
        # Compute interim quantities. 
        ω = ω_weights(grid, baseline.θ, baseline.σ-1)
        value_matching = sol.v_tilde[1] - dot(sol.v_tilde, ω) + sol.x
        free_entry = sol.v_tilde[1] - sol.x*(1-baseline.χ)/baseline.χ
        # Test. 
        @test value_matching < 1e-8
        @test_broken free_entry < 1e-8
    end    
end 