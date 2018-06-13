# Not in the regression test
using PerlaTonettiWaugh, Plots, BenchmarkTools, Sundials, Base.Test

# Testing code
x_min = 0.01
x_max = 1.0
M = 20
T = 10.0
rho = 0.15
sigma_bar = 0.1
c_tilde(t, x) = x + 0.0001*t
sigma_tilde(t, x) =  sigma_bar
mu_tilde(t, x) = 0.1*x *(1.0 + 4.0 * t)
basealgorithm = CVODE_BDF() #CVODE_CDF(linear_solver=:GMRES) #ImplicitEuler() #A reasonable alternative. Algorithms which don't work well: Rosenbrock23(), Rodas4(), KenCarp4()
plotevery = 5

prob = createODEproblem(c_tilde, sigma_tilde, mu_tilde, x_min, x_max, M, T, rho)
sol = solve(prob, basealgorithm)
plot(sol, vars=1:plotevery:M)
@assert(issorted(sol[end]))
# @benchmark solve($prob, $basealgorithm) #Benchmark

#Solve backwards with a DAE and a trivial algebraic equation
#u_ex(1:M) = u(), following the ODE
#u_ex(M+1) = 1.0 #(i.e., a trivial linear function)

probDAE = createDAEproblem(c_tilde, sigma_tilde, mu_tilde, x_min, x_max, M, T, rho)
solDAE = solve(probDAE, IDA())
plot(solDAE, vars=1:plotevery:M)
@show(issorted(solDAE[end][1:M]))

# @benchmark solve($probDAE, IDA())

#Check they are "reasonably" close
@show norm(sol[1] - solDAE[1][1:M])
@show norm(sol[end] - solDAE[end][1:M])
