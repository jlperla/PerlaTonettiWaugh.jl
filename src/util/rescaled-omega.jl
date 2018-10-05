# Given E_hat such that E_hat(T) = 0, perform rescaling to have Ω (and E) such that
# Ω(0) = Ω_0 AND Ω(T) ≈ Ω_T AND E(T) ≈ δ
struct RescaledΩ
    E::Function
    E_hat::Function
    Ω::Function

    function RescaledΩ(E_hat::Function, T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        @assert E_hat(T) ≈ 0
        M = log(Ω_T/Ω_0) / quadgk(E_hat, 0, T)[1]

        Ω_derivative(Ω,p,t) = M*E_hat(t)*Ω
        Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15)
        Ω(t::Float64)::Float64 = t <= T ? Ω_solution(t) : Ω_solution(T)
        E(t::Float64)::Float64 = M*E_hat(t) + δ
        

        new(E, E_hat, Ω)
    end

    # RescaledΩ with piecewise cubic spline interpolation
    function RescaledΩ(E_hat_vec::Array{Float64,1}, T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        @assert length(E_hat_vec) >= 4 # the length of E_hat_vec should be bigger than 3 (cubic spline order)
        @assert E_hat_vec[end] > E_hat_vec[1] # E_hat_vec should have an increasing order, at least for endpoints
        E_hat_vec_range = E_hat_vec[end] - E_hat_vec[1]
        E_hat_vec_scaled = (E_hat_vec .- E_hat_vec[1]) ./ E_hat_vec_range .- 1.0 
        ts = range(0.0, stop=T, length=length(E_hat_vec))
        E_hat_interpolation = Spline1D(ts, E_hat_vec_scaled; k = 3) # cubic spline
        E_hat(t) = E_hat_interpolation(t)

        M = log(Ω_T/Ω_0) / integrate(E_hat_interpolation, 0, T) # (eq:D.2)
        Ω_derivative(Ω,p,t) = M*E_hat(t)*Ω # (eq:16), by replacing E(t) - δ by M*E_hat(t)*Ω from (eq:D.1)
        Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15) # solution for (eq:16)
        Ω(t::Float64)::Float64 = t <= T ? Ω_solution(t) : Ω_solution(T)
        E(t::Float64)::Float64 = M*E_hat(t) + δ # (eq:D.1)

        new(E, E_hat, Ω)
    end
end
(f::RescaledΩ)(t) = f.Ω(t)