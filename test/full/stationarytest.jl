# Define common objects. 
    baselineparams = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 

#=
    Algebraic tests. 
=#

# Run for some different starting poitns (arbitrarily chosen).
    init_x1 = [3.0, -1, -2] # For numerical
    init_x2 = [0.25, 3.0, 1.0] # Breaks around 0.3 for g. 
    init_x3 = [0.0190, 1.434969, 1.06517] # ~ equilibrium value. 

    # Tests for vanilla parameters. 
    # Broken res1_alg = stationary_algebraic(baselineparams(), init_x1)
    res2_alg = stationary_algebraic(baselineparams(), init_x2)
    res3_alg = stationary_algebraic(baselineparams(), init_x3)
    @test res2_alg.g ≈ res3_alg.g atol = atol = 1e-5 

    # Tests for new d. BROKEN FOR NOW. 
    # res1_new = stationary_algebraic(baselineparams(d = 5), init_x1)
    # res2_new = stationary_algebraic(baselineparams(d = 5), init_x2) 
    # # BROKEN res3_new = stationary_algebraic(baselineparams(d = 5), init_x3)
    # @test res1_new.g ≈ res3_new.g atol = 1e-5
    # @test res1_new.g ≈ res2_new.g atol = 1e-5

    # Test for error handling. 
    # @test_throws ErrorException stationary_algebraic(baselineparams(ζ = 0.0001, γ = 0.24))

# Benchmarks.
    # @btime result = stationary_algebraic(x1, baselineparams())
    # @btime result = stationary_algebraic([0.15, 4.0, 2.6], baselineparams()) # Try different params. 
    # Anderson acceleration? Need to find a paramset that works. 

#= 
    Numerical tests. 
=#
    z = unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 250)' range(2.0, stop = 7.0, length = 100)'])
    res1_num = stationary_numerical(baselineparams(), z, init_x1) 
    res2_num = stationary_numerical(baselineparams(), z, init_x2)
    res3_num = stationary_numerical(baselineparams(), z, init_x3)
    res_def_num = stationary_numerical(baselineparams(), z)

    # Residuals
    x = res1_num.x 
    ω = ω_weights(z, baselineparams().θ, baselineparams().σ-1)
    value_matching = v_tilde -> v_tilde[1] - dot(v_tilde, ω) + x
    free_entry = v_tilde -> v_tilde[1] - x*(1-baselineparams().χ)/baselineparams().χ
    @test value_matching(res1_num.v_tilde) ≈ 0.0 atol = 1e-8
    @test free_entry(res1_num.v_tilde) ≈ 0.0 atol = 1e-8

    # Similar to each other. 
    @test res1_num.g ≈ res2_num.g 
    @test res2_num.g ≈ res3_num.g 
    @test res3_num.g ≈ res_def_num.g 

    # Similar-ish to algebraic. 
    @test res1_num.g ≈ res2_alg.g atol = 1e-1