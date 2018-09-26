# make sure that all passes do pass
@testset "PolynomialΩ" begin
    T = 25.0
    Ω_0 = 1.0
    Ω_T = 9.0
    δ = 0.05

    Ω = PolynomialΩ(T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ
    Ω = PolynomialΩ([0.1], T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ
    Ω = PolynomialΩ([0.1; 0.2], T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ
    Ω = PolynomialΩ([-0.1; 0.2; 0.03], T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ
end