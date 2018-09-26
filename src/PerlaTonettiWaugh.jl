module PerlaTonettiWaugh

# Dependencies.
using NLsolve, DifferentialEquations, BandedMatrices, Sundials, Distributions, Roots, Optim, QuantEcon, LinearAlgebra, Random, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, ContinuousTransformations
import Parameters: @with_kw, @unpack

# Utilities
include("util/compactify.jl")
include("util/check-in-interval.jl")
include("util/polynomial-omega.jl")
# General discretization files.
include("diffusionoperators.jl")
include("quadrature.jl")
# Simple model.
include("simple/dynamic.jl")
include("simple/stationary.jl")
include("simple/calculate_residuals.jl")
# Full model.
include("full/stationary.jl") # Includes params.
include("full/dynamic.jl")

export stationary_algebraic, stationary_numerical, simpleODE, simpleDAE, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, calculate_residuals, rescaled_diffusionoperators, diffusionoperators, solve_dynamics, entry_residuals


# export utility for unit tests
export Compactifier, Decompactifier, zero_is_in_interval, is_positive_in_interval, PolynomialΩ

end # module