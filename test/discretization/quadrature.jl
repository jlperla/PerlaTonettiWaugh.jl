using PerlaTonettiWaugh, Distributions


# Test for moments of a truncated normal distribution
testDist1 = TruncatedNormal(0, 1, -1.5, 1.5)
var1 = 0.5515244157615512
pdf1 = x -> pdf(testDist1, x)

z = collect(linspace(-1.5, 1.5, 1000)) # Grid of 1000 uniform points 

omega = irregulartrapezoidweights(z, pdf1)
omega_prime = irregulartrapezoidweights(z, testDist1)
h = x -> (x-0)^2
E_val = dot(h.(z), omega)
E_val_prime = dot(h.(z), omega_prime)
@test E_val ≈ var1 atol = 1e-7
@test E_val_prime ≈ var1 atol = 1e-7

# Test for moments of a truncated exponential distribution
λ = 3
ubound = 3
testDist2 = Truncated(Exponential(λ), 0, ubound)
pdf2 = x -> pdf(testDist2, x)
k = ubound/λ
z = collect(linspace(0, 3.0, 1000))
omega = irregulartrapezoidweights(z, pdf2)
omega_prime = irregulartrapezoidweights(z, testDist2)

    # Test mean
    mean1 = λ * (1 - (k+1)*e^(-k)) / (1 - e^(-k)) # Not ∈ Distributions.jl
    h = x -> x
    E_val = dot(h.(z), omega)
    E_val_prime = dot(h.(z), omega_prime)
    @test E_val ≈ mean1 atol = 1e-6
    @test E_val_prime ≈ mean1 atol = 1e-6

    # Test second raw moment 
    mom2 = 2 * λ^2 * (1 - 0.5*(k^2 + 2*k + 2)*e^(-k))/(1 - exp(-k)) # Not ∈ Distributions.jl
    h = x -> x^2
    E_val = dot(h.(z), omega)
    E_val_prime = dot(h.(z), omega_prime)
    @test E_val ≈ mom2 atol = 1e-6
    @test E_val_prime ≈ mom2 atol = 1e-6

