module PerlaTonettiWaugh

# Dependencies.
using NLsolve, DifferentialEquations, NamedTuples, BandedMatrices, Sundials, Distributions, Roots
import Parameters: @with_kw, @unpack

# General discretization files.
include("diffusionoperators.jl")
include("quadrature.jl")
# Simple model.
include("simple/dynamic.jl")
include("simple/stationary.jl")
include("simple/calculate_residuals.jl")
# Full model.
include("full/stationary.jl") # Includes params.

export f!, simpleODE, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, calculate_residuals, rescaled_diffusionoperators, diffusionoperators

end # module
