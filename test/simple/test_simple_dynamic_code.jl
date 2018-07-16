# a simple test of simple_dynamic_ODE.jl

using PerlaTonettiWaugh, Plots, BenchmarkTools, Sundials, Base.Test,Interpolations, QuantEcon
solve(x,y)=Sundials.solve(x,y)
# Testing code
z_min = 0.01
z_max = 1.0
M = 20
z=linspace(z_min,z_max,M)
T = 10.0
sigma_bar = 0.1
π(t, x) = exp(x)*(1.0+0.0*t)
sigma(t, x) =  sigma_bar
g(t, x) = 0.02 + 0.0*x + 0.0*t
ζ(t, x) = 14.5 + 0.0*x + 0.0*t
r(t, x) = 0.05 + 0.0*x + 0.0*t
gamma=0.005;

params=@NT(gamma=gamma, sigma=sigma, π=π, ζ=ζ, r = r)
settings=@NT(z = z,g = g, T = T)
basealgorithm = CVODE_BDF() #CVODE_CDF(linear_solver=:GMRES) #ImplicitEuler() #A reasonable alternative. Algorithms which don't work well: Rosenbrock23(), Rodas4(), KenCarp4()
plotevery = 5

