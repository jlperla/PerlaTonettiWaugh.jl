module PerlaTonettiWaugh

# package code goes here
using DifferentialEquations, NamedTuples, Parameters, MacroTools, BenchmarkTools, BandedMatrices, Sundials

include("diffusionoperators.jl")
include("simple/ODEalgorithm.jl")
include("simple/DAEalgorithm.jl")
include("simple/algebraicstationary.jl")
include("simple/numericalstationary.jl")
include("full/algebraicstationary.jl")
include("utilities.jl")
export diffusionoperators, simplecreateODEproblem, simplecreateDAEproblem, @with_kw, stationary_algebraic_simple

end # module
