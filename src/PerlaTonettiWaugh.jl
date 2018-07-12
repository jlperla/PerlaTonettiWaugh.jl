module PerlaTonettiWaugh

# package code goes here
using DifferentialEquations, NamedTuples, Parameters, MacroTools, BandedMatrices, Sundials, Distributions

include("diffusionoperators.jl")
include("simple/ODEalgorithm.jl")
include("irregulardiffusionoperators.jl")
include("simple/irregularODEalgo.jl")
include("simple/DAEalgorithm.jl")
include("simple/algebraicstationary.jl")
include("simple/numericalstationary.jl")
include("utilities.jl")
include("full/algebraicstationary.jl")
include("quadrature.jl")
export f!, diffusionoperators, createsimpleODEproblem, createsimpleDAEproblem, @kw_nt, stationary_algebraic_simple, createsimplenonuniformODEproblem, irregulardiffusionoperators, stationary_numerical_simple, irregulartrapezoidweights

end # module
