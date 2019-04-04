module PerlaTonettiWaugh

# Dependencies.
using DifferentialEquations, Sundials, LinearAlgebra, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, QuadGK, NLopt, LeastSquaresOptim, BlackBoxOptim, CSV
using NLSolversBase
using SimpleDifferentialOperators
import Parameters: @with_kw, @unpack

# General utilities files.
include("utils/params.jl")
include("utils/quadrature.jl")
include("utils/solve-system.jl")
include("utils/consumption-equivalent.jl")
include("utils/display_stationary_sol.jl")
# Simple model.
include("simple/dynamic.jl")
include("simple/stationary.jl")
# Full model.
include("full/stationary.jl")
include("full/dynamic.jl")
include("full/transition.jl")

export parameter_defaults, settings_defaults
export params_simple, settings_simple
export solve_simple_dae
export stationary_algebraic, stationary_numerical, simpleDAE, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, solve_dynamics, welfare
export weighted_residuals_given_E_nodes_interior
export solve_full_model
export consumption_equivalent, display_stationary_sol
export f!, f!_simple, f_simple, solve_system

end # module
