# Given E_hat such that E_hat(T) = 0, perform rescaling to have Ω (and E) such that
# Ω(0) = Ω_0 AND Ω(T) ≈ Ω_T AND E(T) ≈ δ
struct RescaledΩ
    E::Function
    E_hat::Function 
    Ω::Function

    # linear case
    function RescaledΩ(E_hat::Function, T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        @assert E_hat(T) ≈ 0
        M = log(Ω_T/Ω_0) / quadgk(E_hat, 0, T)[1]

        E(t) = M*E_hat(t) + δ
        Ω(t) = Ω_0*exp(M*quadgk(E_hat, 0, t)[1]) 

        new(E, E_hat, Ω)
    end
end
(f::RescaledΩ)(t) = f.Ω(t)