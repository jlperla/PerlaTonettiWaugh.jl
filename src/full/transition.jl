#=
  Functions to solve the model (i.e., return equilibrium quantities and a set of E nodes to interpolate) on a transition from one steady state to another.

  The first set of functions are caller functions, which users would interact with. Then come the auxiliary functions used by the callers.
=#

# Global solver.
function solve_full_model_global(settings; impose_E_monotonicity_constraints = true)
  if (settings.transition_iterations < 1)
    E_nodes_interior = settings.transition_x0
    return (solution = solve_model_from_E_nodes(E_nodes_interior, settings; detailed_solution = true), E_nodes = E_nodes_interior)
  end
  ranges = map(i->(settings.transition_lb[i], settings.transition_ub[i]), 1:length(settings.transition_x0))
  ssr(residuals) = sum(residuals.^2)
  result = bboptimize(x -> ssr(weighted_residuals_given_E_nodes_interior(impose_E_monotonicity_constraints ? sort(x) : x, settings));
                      SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = settings.transition_iterations)
  E_nodes_found = best_candidate(result)
  E_nodes_found = impose_E_monotonicity_constraints ? sort(E_nodes_found) : E_nodes_found

  return (solution = solve_model_from_E_nodes(E_nodes_found, settings; detailed_solution = true), E_nodes = E_nodes_found)
end

function solve_full_model(settings; impose_E_monotonicity_constraints = true, write_csv = false, csvpath = nothing)
   # constraint for increasing E nodes
  function constraints_increasing_E!(h, x, jacobian_t)
    L = length(x)
    # A is a matrix whose ith row has 1 in ith col and -1 in (i+1)th col
    # so ith element of A*x imposes x[i] <= x[i+1]
    # note that this imposes x[end] <= 0 as well as the last row is [0; 0; ...; 0; 1].
    A = LinearAlgebra.Tridiagonal(zeros(L-1),ones(L),-ones(L-1))
    if length(jacobian_t) > 0 # transpose of the Jacobian matrix
        jacobian_t[:] = A'
    end
    h[:] = A*x
  end
   result = solve_system(x -> weighted_residuals_given_E_nodes_interior(x, settings), settings.transition_x0;
                    lb = nothing, ub = fill(0.0, length(settings.transition_x0)),
                    constraints_fg! = impose_E_monotonicity_constraints ? constraints_increasing_E! : nothing,
                    algorithm = impose_E_monotonicity_constraints ? :LD_SLSQP : :LD_LBFGS,
                    iterations = settings.transition_iterations)
    solution = solve_model_from_E_nodes(result, settings; detailed_solution = true)

    # output caching
    if write_csv
      CSV.write(csvpath, solution.results)
    end

  return (solution = solution, E_nodes = result)
end

#=
  Auxiliary functions. Organized from most fundamental to least fundamental.
=#

# Returns model solution (i.e., solve_dynamics output) given a set of E nodes. Used by all solvers.
function solve_model_from_E_nodes(E_nodes_interior, settings; detailed_solution = false, interp = LinearInterpolation)
  @unpack T, params_T, stationary_sol_T, Ω_0 = settings
  δ = params_T.δ
  Ω_T = stationary_sol_T.Ω
  # fix the point at T to be zero and sort candidate if needed
  E_nodes = [-1.0; E_nodes_interior; 0.0] # fix the end points to remove indeterminacy problems
  ts = range(0.0, stop=T, length=length(E_nodes))
  E_hat_interpolation = interp(ts, E_nodes) # might worth trying cubic spline
  E_hat(t) = E_hat_interpolation(t)
  E_hat_integral = quadgk(E_hat, 0, T)[1]
  # Formulate and solve ODEProblem
  Q = log(Ω_T/Ω_0) / E_hat_integral # when Ω_T = Ω_0, then Q = 0 so that E(t) is constant with δ as expected
  Ω_derivative(Ω,p,t) = Q*E_hat(t)*Ω
  Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15) # if this fails, error will be thrown
  Ω(t) = t <= T ? Ω_solution(t) : Ω_solution(T)
  E(t) = Q*E_hat(t) + δ
  # solve the dynamics and get the resulting entry_residual vector; if solution is not valid, return Inf
  return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; detailed_solution = detailed_solution)
end

function weighted_residuals_given_E_nodes_interior(E_nodes_interior, settings)
  entry_residuals_nodes_count = length(E_nodes_interior)
  solved = try solve_model_from_E_nodes(E_nodes_interior, settings).results catch; return fill(10e20, entry_residuals_nodes_count) end
  entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
  entry_residuals_nodes = range(0, stop = solved.t[end], length = entry_residuals_nodes_count + 2)

  # return weighted residuals
  return settings.weights .* entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
end
