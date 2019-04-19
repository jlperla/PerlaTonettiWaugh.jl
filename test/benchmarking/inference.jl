# UTILITIES
# Diffusion operators
function testdiffusion()
    ξ = 1.0
    z = 0.:0.01:5.
    bc = (Mixed(ξ), Mixed(ξ))
    a = (L₁₋bc = L₁₋bc(z, bc), L₁₊ = L₁₊bc(z, bc), L₂bc = L₂bc(z, bc))
    b = (L₁₋bc = L₁₋bc(collect(z), bc), L₁₊ = L₁₊bc(collect(z), bc), L₂bc = L₂bc(collect(z), bc)) # test both the AbstractRange and StepRangeLen methods
    return (a = a, b = b)
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
    params = parameters_simple() # Arbitrary parameter set
    settings = settings_simple()
    PerlaTonettiWaugh.simpleDAE(params, settings) # non-inferrable because of the constructor call
end

println("Benchmarking simpleDAE constructor")
@btime test_simpleDAE()
