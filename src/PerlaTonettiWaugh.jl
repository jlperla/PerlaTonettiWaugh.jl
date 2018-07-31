module PerlaTonettiWaugh

# Dependencies. 
using NLsolve, DifferentialEquations, NamedTuples, BandedMatrices, Sundials, Distributions, Roots
import Parameters: @with_kw, @unpack

# General discretization files. 
include("diffusionoperators.jl")
include("quadrature.jl")
# Simple model.
include("simple/diffeqs.jl")
include("simple/algebraicstationary.jl")
include("simple/numericalstationary.jl")
include("simple/calculate_residuals.jl")
# Full model. 
include("full/algebraicstationary.jl") # Includes params.

export f!, diffusionoperators, simpleODEproblem, simpleDAEproblem, simpledynamicODEproblem, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, irregulartrapezoidweights, calculate_residuals, rescaled_diffusionoperators

end # module
