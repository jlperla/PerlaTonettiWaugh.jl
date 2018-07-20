# This is a test of calculate_residual.jl

using PerlaTonettiWaugh, BenchmarkTools, Sundials, Base.Test,Interpolations, QuantEcon, NamedTuples

# setting up
z_min = 0.0
z_max = 5.0
M = 100
N = 50
z=linspace(z_min,z_max,M)

T = 10.0
t = linspace(0, T, N)

π(t, x) = exp(x)*(1.0+0.0*t);
σ_n = 0.02
α = 2.1
ζ_n = 14.5
r_n = 0.05
ζ(t) = ζ_n + 0.0*t
r(t, x) = r_n + 0.0*x + 0.0*t
σ(t, x) = σ_n + 0.0*x + 0.0*t
γ=0.005;

params=@NT(γ=γ, σ=σ, π=π, ζ=ζ, r=r, α=α)
params_n=@kw_nt(γ=0.005, σ=0.02, α=2.1, r=0.05, ζ=14.5)

# solve for numerical g_T as test
res=stationary_numerical_simple(params_n(),z)
g_analytic=res.g
v_analytic=res.v
@show v_analytic

# Create the interpolation object of g
g_vector = g_analytic + 0.01 * t
g_int=LinInterp(t, g_vector)
g(t, x) = g_int(t) + 0.0 * x

resid =calculate_residuals(params, g, z, T)
@show resid