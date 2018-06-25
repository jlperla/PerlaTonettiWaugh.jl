#Just implements the analytical solution we have in the notes returning back functions and g

### Start with substitutions from Appendix H 
F(z::Real) = 1 - z^(-θ); # (H.1)
S = θ(g - μ - θ*v^2*0.5); # (H.2)
nu = (μ-g)/v^2 + √(((g-μ)/v^2)^2 + (r-g)/(0.5 * v^2)); # (H.3)
a = (r - g - (σ - 1)*(μ - g + (σ -1)*v^2/2)); # (H.4)
b = (1 - a(r-g))*d^(1 - σ)̂*̂z^(nu + σ - 1); # (H.5)
r = ρ + γg + δ; # (H.6)
L̃ = Ω((N-1)(1 - F(ẑ))κ + (1-η)ζ(S+δ/χ)); # (H.7)
# (H.8)
̂z = d(κ/̄πₘ)^(1/(σ - 1)) # (H.9)
w = (1/̄σ)*̄z; # (H.10)
x = ζ(1-η + ηΘ/w); # (H.11)