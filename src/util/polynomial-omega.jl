# polynomials such that
# Ω(0) = Ω_0 AND Ω(T) ≈ Ω_T AND E(T) ≈ δ
struct PolynomialΩ
    E::Function
    E_derivative::Function
    Ω::Function

    # linear case
    function PolynomialΩ(T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        Ω = t -> Ω_0 * (Ω_0/Ω_T)^(t*(t-2.0*T)/(T*T))
        E = t -> 2*(t-T)*log(Ω_0/Ω_T)/T^2 + δ
        E_derivative = t -> 2*(t-T)*log(Ω_0/Ω_T)/T^2 + δ

        @assert Ω(0) ≈ Ω_0
        @assert Ω(T) ≈ Ω_T
        @assert E(T) ≈ δ

        new(E, E_derivative, Ω)
    end

    function PolynomialΩ(c::Array{Float64,1}, T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        deg = length(c) + 1 
        if (deg < 1 || deg > 5) throw("Up to quartic polynomials are supported.") end
        # note that deg == 1 implies that c is empty (which will return type error as c is a Float64 array)
        if (deg == 2)
            Ω = t -> exp(c[1]*t*((t-T)^2)/3 ) * Ω_0 * (Ω_0/Ω_T)^(t*(t-2.0*T)/(T*T))
            E = t -> c[1]*(3*t*t-4*t*T+T*T)/3 + 2*(t-T)*log(Ω_0/Ω_T)/T^2 + δ
            E_derivative = t -> c[1]*(6*t-4*T)/3 + 2*log(Ω_0/Ω_T)/T^2
        end
        if (deg == 3)
            Ω = t -> exp( t*((t-T)^2)*(4*c[1]+3*c[2]*(t+2*T)) /12 ) * Ω_0 * (Ω_0/Ω_T)^(t*(t-2.0*T)/(T*T))
            E = t -> (12*t*log(Ω_0/Ω_T) + T*( (t-T)*T*(6*t*(c[1]+c[2]*t) - 2*(c[1]-3*c[2]*t)*T - 3*c[2]*T*T)  + 12*log(Ω_T/Ω_0) ) ) / (6*T*T) + δ
            E_derivative = t -> 2*c[1]*t+3*c[2]*t^2-4*c[1]*T/3-3*c[2]*T^2/2+2*log(Ω_0/Ω_T)/T^2
        end 
        if (deg == 4)
            Ω = t -> exp( t*(t-T)^2*(20*c[1]+15*c[2]*(t+2T)+12*c[3]*(t^2+2*t*T+3T^2)) ) * Ω_0 * (Ω_0/Ω_T)^(t*(t-2.0*T)/(T*T) )
            E = t -> t^2*(c[1]+t*(c[2]+c[3]*t))-(4*c[1]*t*T)/3+(2*c[1]-9*c[2]*t)T^2/6+(5*c[2]-16*c[3]*t)*T^3/10+(3*c[3]*T^4)/5+δ+(2(t-T)*log(Ω_0/Ω_T))/T^2
            E_derivative = t -> 2*c[1]*t + 3*c[2]*t^2 + 4*c[3]*t^3 - (4*c[1]*T)/3 - (3*c[2]*T^2)/2 - (8*c[3]*T^3)/5 + (2*log(Ω_0/Ω_T))/T^2
        end

        @assert Ω(0) ≈ Ω_0
        @assert Ω(T) ≈ Ω_T
        @assert E(T) ≈ δ

        new(E, E_derivative, Ω)
    end
end

(f::PolynomialΩ)(t) = f.Ω(t)