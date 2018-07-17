# a simple test of simple_dynamic_ODE.jl

using PerlaTonettiWaugh, Plots, BenchmarkTools, Sundials, Base.Test,Interpolations, QuantEcon, NamedTuples
solve(x,y)=Sundials.solve(x,y)
# setting up
z_min = 0.01
z_max = 5.0
M = 100
z=linspace(z_min,z_max,M)
T = 10.0
π(t, x) = exp(x)*(1.0+0.0*t);
σ_n = 0.02
α = 2.1
ζ_n = 14.5
r_n = 0.05
ζ(t, x) = ζ_n + 0.0*x + 0.0*t
r(t, x) = r_n + 0.0*x + 0.0*t
σ(t, x) = σ_n + 0.0*x + 0.0*t
γ=0.005;

# solve for numerical g_T as test
params_n=@kw_nt(γ=0.005, σ=0.02, α=2.1, r=0.05, ζ=14.5)

res=stationary_numerical_simple(params_n(),z)
g_analytic=res.g
v_analytic=res.v
@show v_analytic

g(t, x) = g_analytic + 0.0*x + 0.0*t;

params=@NT(γ=γ, σ=σ, π=π, ζ=ζ, r = r, α=α)
settings=@NT(z = z,g = g, T = T)
basealgorithm = CVODE_BDF() #CVODE_CDF(linear_solver=:GMRES) #ImplicitEuler() #A reasonable alternative. Algorithms which don't work well: Rosenbrock23(), Rodas4(), KenCarp4()
plotevery = 5

# 1. test starting from g_analytic, v_T are close to analytical result
prob = create_dynamic_ODE(params,settings)
sol = solve(prob, basealgorithm)
# plot(sol, vars=1:plotevery:M)
@test issorted(sol[end])
@show sol[end]


# 2. test whether ODE result close to analytic
tol=1e-8
@test norm(v_analytic-sol[end])<tol

# If start at some time varying function of g(t) (change g grid), notice we shouldn't expect this is close to v_analytic

g_g(t, x) = g_analytic + 0.0*x + 0.01*t;

settings_g=@NT(z = z,g = g_g, T = T)

prob_g = create_dynamic_ODE(params,settings_g)
sol_g = solve(prob_g, basealgorithm)
@show sol_g[end]
@test issorted(sol_g[end])

# 3. Test non uniform grid z

z_add=linspace(z_min,z_max,37)
z_comb=unique(sort([z;z_add]))

settings_ir=@NT(z = z_comb,g = g, T = T)

prob_ir = create_dynamic_ODE(params,settings_ir)
sol_ir = solve(prob_ir, basealgorithm)
@show sol_ir[end]
@test issorted(sol_ir[end])

# interpolate uniform onto non uniform grid
sol_int=LinInterp(z, sol[end])
@show norm(sol_int.(z_comb)-sol_ir[end])
@show norm(sol_ir[end,end]-sol[end,end])<1e-2
