using PerlaTonettiWaugh, Test, LinearAlgebra, Statistics, Compat
using Distributions, Sundials, BenchmarkTools, Parameters, QuantEcon, Interpolations, NLsolve, Optim, DifferentialEquations, DiffEqCallbacks, Random, Statistics, Interpolations, Roots, QuadGK

@elapsed begin 
  @time @testset "Simple Stationary" begin include("simple/stationarytest.jl") end
  @time @testset "Full Stationary" begin include("full/stationarytest.jl") end
  @time @testset "Residuals/Dynamic ODE Tests" begin include("simple/residualstest.jl") end
  @time @testset "Full Dynamic" begin include("full/dynamictest.jl") end
  @time @testset "Entry residuals" begin include("full/entryresidualstest.jl") end
end 
