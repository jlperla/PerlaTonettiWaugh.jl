# make sure that all passes do pass
@testset "RescaledΩ" begin
    T = 25.0
    Ω_0 = 1.0
    Ω_T = 9.0
    δ = 0.05

    # linear E_hat
    Ω = RescaledΩ(t -> (t-T), T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ

    # nonlinear E_hat
    Ω = RescaledΩ(t -> -t*(t-T), T, Ω_0, Ω_T, δ)
    @assert Ω(0.0) ≈ Ω_0
    @assert Ω(T) ≈ Ω_T
    @assert Ω.E(T) ≈ δ
end