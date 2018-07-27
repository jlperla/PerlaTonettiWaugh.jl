using PerlaTonettiWaugh, Base.Test
using Distributions, Sundials, BenchmarkTools, QuantEcon, Interpolations, Parameters, NamedTuples, NLsolve

tic()
@time @testset "Analytical Simple Stationary" begin include("simple/algebraicstationarytest.jl") end
@time @testset "Numerical Simple Stationary" begin include("simple/numericalstationarytest.jl") end
@time @testset "Analytical Full Stationary" begin include("full/algebraicstationarytest.jl") end
@time @testset "Discretization Tests" begin include("discretization/discretizationtest.jl") end 
@time @testset "Quadrature Tests" begin include("discretization/quadrature.jl") end 
@time @testset "Residuals/Dynamic ODE Tests" begin include("simple/residualstest.jl") end 
@time @testset "ODE and DAE Tests" begin include("simple/diffeqtest.jl") end
toc()



