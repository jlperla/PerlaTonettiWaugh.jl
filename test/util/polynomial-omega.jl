# make sure that all passes do pass
@testset "PolynomialΩ" begin
    T = 25.0
    Ω_0 = 1.0
    Ω_T = 9.0
    δ = 0.05
    p = PolynomialΩ([0.1], T, Ω_0, Ω_T, δ)
    @assert p.Ω(0.0) ≈ Ω_0
    @assert p.Ω(T) ≈ Ω_T
    @assert p.E(T) ≈ δ
    p = PolynomialΩ([0.1; 0.2], T, Ω_0, Ω_T, δ)
    @assert p.Ω(0.0) ≈ Ω_0
    @assert p.Ω(T) ≈ Ω_T
    @assert p.E(T) ≈ δ
    p = PolynomialΩ([-0.1; 0.2; 0.03], T, Ω_0, Ω_T, δ)
    @assert p.Ω(0.0) ≈ Ω_0
    @assert p.Ω(T) ≈ Ω_T
    @assert p.E(T) ≈ δ
end