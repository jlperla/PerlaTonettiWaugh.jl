using PerlaTonettiWaugh, Distributions

# Test for moments of a truncated normal distribution
testDist1 = TruncatedNormal(0, 1, -1.5, 1.5)
var1 = 0.5515244157615512
pdf1 = x -> pdf(testDist1, x)

z = collect(linspace(-1.5, 1.5, 1000)) # Grid of 1000 uniform points 

omega = irregulartrapezoidweights(z, pdf1)
h = x -> (x-0)^2
E_val = dot(h.(z), omega)
@test E_val ≈ var1 atol = 1e-7

# Test for moments of a truncated exponential distribution
λ = 3
ubound = 3
testDist2 = Exponential(λ)
pdf2 = x -> x <= ubound && x > 0 ? (inv(λ) * exp(-x * inv(λ)))/(1 - exp(-1 * ubound/λ)) : 0.0
k = ubound/λ
mean1 = λ * (1 - (k+1)*e^(-k)) / (1 - e^(-k))
z = collect(linspace(0, 3.0, 1000))
omega = irregulartrapezoidweights(z, pdf2)

    # Test mean
    mean1 = 3 * (1 - 2e^(-1)) / (1 - e^(-1))
    h = x -> x
    E_val = dot(h.(z), omega)
    @test E_val ≈ mean1 atol = 1e-6

    # Test second raw moment 
    mom2 = 2 * λ^2 * (1 - 0.5*(k^2 + 2*k + 2)*e^(-k))/(1 - exp(-k))
    h = x -> x^2
    E_val = dot(h.(z), omega)
    @test E_val ≈ mom2 atol = 1e-6

