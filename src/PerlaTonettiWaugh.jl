module PerlaTonettiWaugh

# Dependencies.
using DifferentialEquations, Sundials, LinearAlgebra, DataFrames, DataFramesMeta, DiffEqCallbacks, Interpolations, QuadGK, CSV
using NLsolve
using SimpleDifferentialOperators
using JSON
import Parameters: @with_kw, @unpack

# General utilities files.
include("utils/quadrature.jl")
# Simple model.
include("simple/params.jl")
include("simple/dynamic.jl")
include("simple/stationary.jl")
# Full model.
include("full/params.jl")
include("full/static.jl")
include("full/stationary.jl")
include("full/dynamic.jl")

export parameter_defaults, settings_defaults
export parameters_simple, settings_simple
export solve_simple_transition
export stationary_algebraic, stationary_numerical, stationary_algebraic_simple, stationary_numerical_simple, ω_weights, solve_dynamics, welfare
export compare_steady_states
export solve_transition, prepare_results
export consumption_equivalent
export f!, f!_simple

end # module
