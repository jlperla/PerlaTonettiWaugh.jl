using PerlaTonettiWaugh, Base.Test

tic()
@time @testset "Analytical Simple Steady State" begin include("simple/algebraicstationarytest.jl") end
@time @testset "Numerical Simple Steady State" begin include("simple/numericalsteadystatetest.jl") end
@time @testset "Regression Simple Steady State" begin include("simple/regressiontest.jl") end 
@time @testset "Analytical Full Steady State" begin include("full/algebraicstationarytest.jl") end
@time @testset "Discretization Tests" begin include("discretization/discretizationtest.jl") end 
toc()