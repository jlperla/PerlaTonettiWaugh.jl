module PerlaTonettiWaugh

# package code goes here
using DifferentialEquations, NamedTuples, Parameters, MacroTools, BandedMatrices, Sundials

include("diffusionoperators.jl")
include("simple/ODEalgorithm.jl")
include("irregulardiffusionoperators.jl")
include("simple/irregularODEalgo.jl")
include("simple/DAEalgorithm.jl")
include("simple/algebraicstationary.jl")
include("simple/numericalstationary.jl")
include("full/algebraicstationary.jl")
include("utilities.jl")

export diffusionoperators, createsimpleODEproblem, createsimpleDAEproblem, @with_kw, stationary_algebraic_simple, createsimplenonuniformODEproblem, irregulardiffusionoperators

end # module
