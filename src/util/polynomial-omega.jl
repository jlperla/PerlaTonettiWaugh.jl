# polynomials such that
# Ω(0) = Ω_0 AND Ω(T) ≈ Ω_T AND E(T) ≈ δ
struct PolynomialΩ
    E::Function
    E_derivative::Function
    Ω::Function

    function PolynomialΩ(c_and_T::Array{Float64,1}, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        if (length(c_and_T) == 1)
            PolynomialΩ(c_and_T[1], Ω_0, Ω_T, δ)
        else
            PolynomialΩ(c_and_T[1:(end-1)], c_and_T[end], Ω_0, Ω_T, δ)
        end
    end

    # linear case
    function PolynomialΩ(T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        Ω = t -> Ω_0 * (Ω_0/Ω_T)^(t*(t-2.0*T)/(T*T))
        E = t -> 2*(t-T)*log(Ω_0/Ω_T)/T^2 + δ
        E_derivative = t -> 2*log(Ω_0/Ω_T)/T^2

        @assert Ω(0) ≈ Ω_0
        @assert Ω(T) ≈ Ω_T
        @assert E(T) ≈ δ

        new(E, E_derivative, Ω)
    end

    function PolynomialΩ(c::Array{Float64,1}, T::Float64, Ω_0::Float64, Ω_T::Float64, δ::Float64)
        deg = length(c) + 1 
        if (deg < 1 || deg > 6) throw("Up to sextic polynomials are supported.") end
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
        if (deg == 5) # quintic
            Ω = t -> exp((t*(60t*log(Ω_0) - 60t*log(Ω_T) + T*((t - T)^2*T*(20*c[1] + 15*c[2]*t + 12*c[3]*t^2 + 10*c[4]*t^3 + 30*c[2]*T + 24*c[3]*t*T + 20*c[4]*t^2*T + 36*c[3]*T^2 + 30*c[4]*t*T^2 + 40*c[4]*T^3) + 120*log(Ω_T/Ω_0))))/(60*T^2))*Ω_0
            E = t -> c[1]*t^2 + c[2]*t^3 + c[3]*t^4 + c[4]*t^5 + (c[1]*T^2)/3 + (c[2]*T^3)/2 + (3*c[3]*T^4)/5 + (2*c[4]*T^5)/3 + δ - (t*(T^3*(40*c[1] + 45*c[2]*T + 48*c[3]*T^2 + 50*c[4]*T^3) - 60*log(Ω_0) + 60*log(Ω_T)))/(30*T^2) + (2*log(Ω_T/Ω_0))/T
            E_derivative = t -> 2*c[1]*t + 3*c[2]*t^2 + 4*c[3]*t^3 + 5*c[4]*t^4 - (4*c[1]*T)/3 - (3*c[2]*T^2)/2 - (8*c[3]*T^3)/5 - (5*c[4]*T^4)/3 + (2*log(Ω_0/Ω_T))/T^2
        end
        if (deg == 6) # sextic
            Ω = t -> exp((t*(420t*log(Ω_0/Ω_T)+T*((t-T)^2*T*(140*c[1]+105*c[2]*(t+2T)+84*c[3]*(t^2+2t*T+3*T^2)+70*c[4]*(t^3+2*t^2*T+3*t*T^2+4*T^3)+60*c[5]*(t^4+2*t^3*T+3*t^2*T^2+4*t*T^3+5*T^4))+840*log(Ω_T/Ω_0))))/(420*T^2))*Ω_0
            E = t -> c[1]*t^2 + c[2]*t^3 + c[3]*t^4 + c[4]*t^5 + c[5]*t^6 + (c[1]*T^2)/3 + (c[2]*T^3)/2 + (3*c[3]*T^4)/5 + (2*c[4]*T^5)/3 + (5*c[5]*T^6)/7 - 1/210*t*T*(280*c[1] + T*(315*c[2] + 336*c[3]*T + 350*c[4]*T^2 + 360*c[5]*T^3)) + δ + (2*t*log(Ω_0/Ω_T))/T^2 + (2*log(Ω_T/Ω_0))/T
            E_derivative = t -> 2*c[1]*t + 3*c[2]*t^2 + 4*c[3]*t^3 + 5*c[4]*t^4 + 6*c[5]*t^5 - (4*c[1]*T)/3 - (3*c[2]*T^2)/2 - (8*c[3]*T^3)/5 - (5*c[4]*T^4)/3 - (12*c[5]*T^5)/7 + (2*log(Ω_0/Ω_T))/T^2
        end

        @assert Ω(0) ≈ Ω_0
        @assert Ω(T) ≈ Ω_T
        @assert E(T) ≈ δ

        new(E, E_derivative, Ω)
    end
end
(f::PolynomialΩ)(t) = f.Ω(t)