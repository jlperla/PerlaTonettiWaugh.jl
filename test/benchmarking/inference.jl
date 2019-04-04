# UTILITIES
# Diffusion operators
function testdiffusion()
    ξ = 1.0
    z = 0.:0.01:5.
    bc = (Mixed(ξ), Mixed(ξ))
    (L₁₋bc = L₁₋bc(z, bc), L₁₊ = L₁₊(z, bc), L₂bc = L₂bc(z, bc)),
    (L₁₋bc = L₁₋bc(collect(z), bc), L₁₊ = L₁₊(collect(z), bc), L₂bc = L₂bc(collect(z), bc)) # test both the AbstractRange and StepRangeLen methods
end

@test @inferred testdiffusion() == testdiffusion()

# ω_weights
function test_ω()
    z = 0.:0.01:5.
    α = 2.1
    ξ = 1.
    ω_weights(z, α, ξ), ω_weights(collect(z), α, ξ), ω_weights(z, α, 1), ω_weights(collect(z), α, 1), ω_weights(z, 1.0, 2.0) # Different type signatures/values for input
end

@test @inferred test_ω() == test_ω()

# MODEL OBJECTS
# simpleDAE
function test_simpleDAE()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z_ex = 0.:0.1:3., T = 10., g = x -> 0.02)
    simpleDAE(params, settings) # non-inferrable because of the constructor call
end

println("Benchmarking simpleDAE constructor")
@btime test_simpleDAE()
