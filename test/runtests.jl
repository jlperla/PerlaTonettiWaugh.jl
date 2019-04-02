using PerlaTonettiWaugh
using Test, LinearAlgebra, Statistics, Compat
using Sundials, BenchmarkTools, Interpolations, Parameters, SimpleDifferentialOperators

@elapsed begin
  # @time @testset "Simple Stationary" begin include("simple/stationarytest.jl") end
  @time @testset "Full Stationary" begin include("full/stationarytest.jl") end
  # @time @testset "Residuals/Dynamic ODE Tests" begin include("simple/residualstest.jl") end
  @time @testset "Full Dynamic" begin include("full/dynamictest.jl") end
  @time @testset "Entry residuals" begin include("full/entryresidualstest.jl") end
  # @time @testset "Type Stability and Speed" begin include("benchmarking/inference.jl") end
  @time @testset "Utilities" begin include("utils/solve-system.jl") end
end
