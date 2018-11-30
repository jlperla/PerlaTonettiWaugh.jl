# UTILITIES
# Diffusion operators
function testdiffusion()
    ξ = 1.0
    z = 0.:0.01:5.
    rescaled_diffusionoperators(z, ξ), rescaled_diffusionoperators(collect(z), ξ) # test both the AbstractRange and StepRangeLen methods
end

@inferred testdiffusion()

# ω_weights
function test_ω()
    z = 0.:0.01:5.
    α = 2.1
    ξ = 1.
    ω_weights(z, α, ξ), ω_weights(collect(z), α, ξ), ω_weights(z, α, 1), ω_weights(collect(z), α, 1) # Different type signatures for input
end

@inferred test_ω()

# DYNAMIC OBJECTS
# simpleODE
function test_simpleODE()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10., g = x -> 0.02)
    @btime simpleODE($params, $settings)
end

test_simpleODE()

# simpleDAE
function test_simpleDAE()
    params = (μ = 0.0, υ = 0.1, θ = 2.1, r = x -> 1., x = n -> 1., ξ = 1., π_tilde = (t, z) -> 1.) # Arbitrary parameter set
    settings = (z = 0.:0.1:3., T = 10., g = x -> 0.02)
    @btime simpleDAE($params, $settings)
end

test_simpleDAE()
