module PerlaTonettiWaugh

# Dependencies.
using DifferentialEquations, Sundials, Distributions, Roots, QuantEcon, LinearAlgebra, Random, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, QuadGK, NLopt, ForwardDiff, SoftGlobalScope
using LeastSquaresOptim
using BlackBoxOptim
import Parameters: @with_kw, @unpack
import DFOLS

# General utilities files.
include("params.jl")
include("diffusionoperators.jl")
include("quadrature.jl")
include("utils/find-zero.jl")
include("utils/consumption-equivalent.jl")
include("utils/display_stationary_sol.jl")
# Simple model.
include("simple/dynamic.jl")
include("simple/stationary.jl")
include("simple/calculate_residuals.jl")
# Full model.
include("full/stationary.jl") # Includes params.
include("full/dynamic.jl")
include("full/solve-full-model.jl")

export parameter_defaults, settings_defaults
export stationary_algebraic, stationary_numerical, simpleODE, simpleDAE, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, calculate_residuals, rescaled_diffusionoperators, diffusionoperators, solve_dynamics, entry_residuals, welfare
export minimize_residuals, minimize_residuals_python
export solve_full_model_global, solve_full_model_python, solve_continuation
export consumption_equivalent, display_stationary_sol
export f!, f!_simple, f_simple

end # module
