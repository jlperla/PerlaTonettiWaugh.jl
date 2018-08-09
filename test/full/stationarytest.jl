# Define common objects. 
    baselineparams = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 

# Run for some different starting poitns (arbitrarily chosen).
    x1 = [0.2, 2.3, 0.4]
    x2 = [0.25, 3.0, 1.0] # Breaks around 0.3 for g. 
    x3 = [0.0260, 3.7504, 1.0615] # Old equilibrium value. 

    # Tests for vanilla parameters. 
    res1_vanilla = stationary_algebraic_full(x1, baselineparams())
    res2_vanilla = stationary_algebraic_full(x2, baselineparams())
    res3_vanilla = stationary_algebraic_full(x3, baselineparams())
    @test res1_vanilla.g ≈ res2_vanilla.g atol = 1e-5
    @test res2_vanilla.g ≈ res3_vanilla.g atol = atol = 1e-5 

    # Tests for new d. 
    res1_new = stationary_algebraic_full(x1, baselineparams(d = 5))
    @test_throws AssertionError stationary_algebraic_full(x2, baselineparams(d = 5)) # Solver goes into the invalid param region. 
    res3_new = stationary_algebraic_full(x3, baselineparams(d = 5))
    @test res1_new.g ≈ res3_new.g atol = 1e-5


# Benchmarks.
    # @btime result = stationary_algebraic_full(x1, baselineparams())
    # @btime result = stationary_algebraic_full([0.15, 4.0, 2.6], baselineparams()) # Try different params. 
    # Anderson acceleration? Need to find a paramset that works. 