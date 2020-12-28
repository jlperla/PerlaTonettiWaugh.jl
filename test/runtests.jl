using PerlaTonettiWaugh
using Test, LinearAlgebra, Statistics, Compat
using Sundials, BenchmarkTools, Interpolations, Parameters, SimpleDifferentialOperators
using DataFrames, DataFramesMeta

@elapsed begin
  @time @testset "Simple Stationary" begin include("simple/stationarytest.jl") end
  @time @testset "Full Stationary" begin include("full/stationarytest.jl") end
  @time @testset "Simple Dynamic" begin include("simple/dynamictest.jl") end
  @time @testset "Full Dynamic" begin include("full/dynamictest.jl") end
  @time @testset "Type Stability and Speed" begin include("benchmarking/inference.jl") end
end
