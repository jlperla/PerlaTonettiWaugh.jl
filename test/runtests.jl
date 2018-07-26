using PerlaTonettiWaugh, Base.Test
using Distributions, Sundials, BenchmarkTools, QuantEcon, Interpolations, Parameters, NamedTuples

tic()
#@time @testset "Analytical Simple Stationary" begin include("simple/algebraicstationarytest.jl") end
#@time @testset "Numerical Simple Stationary" begin include("simple/numericalstationarytest.jl") end
#@time @testset "Regression Simple Stationary" begin include("simple/regressiontest.jl") end 
#@time @testset "Analytical Full Stationary" begin include("full/algebraicstationarytest.jl") end
#@time @testset "Discretization Tests" begin include("discretization/discretizationtest.jl") end 
@time @testset "ODE and DAE Tests" begin include("simple/diffeqtest.jl") end
#@time @testset "Quadrature Tests" begin include("discretization/quadrature.jl") end 
@time @testset "Residuals Tests" begin include("simple/residualstest.jl") end
toc()


