module PerlaTonettiWaugh

# package code goes here
using DifferentialEquations, NamedTuples, BenchmarkTools, BandedMatrices, Sundials

include("diffusionoperators.jl")
include("simple/ODEalgorithm.jl")
include("simple/DAEalgorithm.jl")
include("simple/analyticalstationary.jl")
include("simple/numericalstationary.jl")
include("full/analyticalstationary.jl")
export diffusionoperators, simplecreateODEproblem, simplecreateDAEproblem

end # module
