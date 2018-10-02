# make sure that all passes do pass
@testset "RescaledΩ" begin
    T = 25.0
    Ω_0 = 1.0
    Ω_T = 9.0
    δ = 0.05

    # linear E_hat
    Ω = RescaledΩ(t -> (t-T), T, Ω_0, Ω_T, δ)
    @test Ω(0.0) ≈ Ω_0
    @test Ω(T) ≈ Ω_T atol = 1e-8
    @test Ω.E(T) ≈ δ

    # nonlinear E_hat
    Ω = RescaledΩ(t -> -t*(t-T), T, Ω_0, Ω_T, δ)
    @test Ω(0.0) ≈ Ω_0
    @test Ω(T) ≈ Ω_T atol = 1e-8
    @test Ω.E(T) ≈ δ
end

@testset "RescaledΩ, convenience constructors" begin
    T = 25.0
    Ω_0 = 1.0
    Ω_T = 9.0
    δ = 0.05

    # RescaledΩ with piecewise linear interpolation E
    Ω = RescaledΩ([-1.0, 0.5, 0.1], T, Ω_0, Ω_T, δ)
    @test Ω(0.0) ≈ Ω_0
    @test Ω(T) ≈ Ω_T atol = 1e-4
    @test Ω.E(T) ≈ δ

    # RescaledΩ with piecewise linear interpolation E
    Ω = RescaledΩ([-1.0, 0.5, 0.1, -0.1], T, Ω_0, Ω_T, δ)
    @test Ω(0.0) ≈ Ω_0
    @test Ω(T) ≈ Ω_T atol = 1e-4
    @test Ω.E(T) ≈ δ
end