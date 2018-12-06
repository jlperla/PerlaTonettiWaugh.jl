module PerlaTonettiWaugh

# Dependencies.
using NLsolve, DifferentialEquations, BandedMatrices, Sundials, Distributions, Roots, QuantEcon, LinearAlgebra, Random, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, QuadGK, Dierckx, NLopt, ForwardDiff
using LeastSquaresOptim
using BlackBoxOptim
import Parameters: @with_kw, @unpack
import DFOLS

# General utilities files.
include("params.jl")
include("diffusionoperators.jl")
include("quadrature.jl")
# Simple model.
include("simple/dynamic.jl")
include("simple/stationary.jl")
include("simple/calculate_residuals.jl")
# Full model.
include("full/stationary.jl") # Includes params.
include("full/dynamic.jl")
include("full/solve-full-model.jl")

export parameter_defaults
export stationary_algebraic, stationary_numerical, simpleODE, simpleDAE, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, calculate_residuals, rescaled_diffusionoperators, diffusionoperators, solve_dynamics, entry_residuals, welfare
export minimize_residuals, minimize_residuals_python
export solve_full_model_global

end # module
