using PerlaTonettiWaugh, Test, LinearAlgebra, Statistics, Compat
using Distributions, Sundials, BenchmarkTools, Parameters, QuantEcon, Interpolations, Optim, DifferentialEquations, DiffEqCallbacks, Random, Statistics, Interpolations, QuadGK

@elapsed begin
  @time @testset "Simple Stationary" begin include("simple/stationarytest.jl") end
  @time @testset "Full Stationary" begin include("full/stationarytest.jl") end
  @time @testset "Residuals/Dynamic ODE Tests" begin include("simple/residualstest.jl") end
  @time @testset "Full Dynamic" begin include("full/dynamictest.jl") end
  @time @testset "Entry residuals" begin include("full/entryresidualstest.jl") end
  @time @testset "Type Stability and Speed" begin include("benchmarking/inference.jl") end
  @time @testset "Utils" begin include("utils/find-zero.jl") end
end
