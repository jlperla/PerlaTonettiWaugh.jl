module PerlaTonettiWaugh

# Dependencies.
using NLsolve, DifferentialEquations, BandedMatrices, Sundials, Distributions, Roots, QuantEcon, LinearAlgebra, Random, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, QuadGK, Dierckx, NLopt, ForwardDiff
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
include("full/dynamic.jl")

export stationary_algebraic, stationary_numerical, simpleODE, simpleDAE, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, calculate_residuals, rescaled_diffusionoperators, diffusionoperators, solve_dynamics, entry_residuals, welfare

end # module
