#=
  Functions to solve the model (i.e., return equilibrium quantities and a set of E nodes to interpolate) on a transition from one steady state to another.

  The first set of functions are caller functions, which users would interact with. Then come the auxiliary functions used by the callers.
=#

# Continuation solver. Essentially traces out a smooth path from some d_T to d_0, using the result from one step as the x0 for the next.
# With default parameters we get: [-1.1965566608914664,-0.7506441878495663,-0.6160629809004126,-0.4490604063250066,-0.3457782022324279,-0.26611899012195817,-0.1703597837124977,-0.12538756277279978,-0.10765936098467942,-0.057163395444530966,-0.05297636889833756,-0.032260206198113685,-0.028143386570548785,-0.04434566417370368]

#
function solve_continuation(d_0, d_T; step = 0.005, params = parameter_defaults(), settings = settings_defaults(), solver = solve_full_model_dfols)
  params_0 = merge(params, (d = d_T,)) # parameters to be used at t = 0
  params_T = merge(params, (d = d_T,)) # parameters to be used at t = T
  z_grid = settings.z
  stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
  stationary_sol_T = stationary_numerical(params_T, z_grid) # solution at t = T
  Ω_0 = stationary_sol_0.Ω;
  settings = merge(settings, (params_T = params_T, stationary_sol_T = stationary_sol_T, Ω_0 = Ω_0, transition_x0 = settings.continuation_x0));
  tempd_0 = params_0.d
  result = 0 # to be overwritten by loop
   while tempd_0 <= d_0
    tempd_0 += step
    params_0 = merge(params, (d = tempd_0,))
    Ω_0 = stationary_numerical(params_0, z_grid).Ω
    settings = merge(settings, (Ω_0 = Ω_0,))
    result = solver(settings)
    settings = merge(settings, (transition_x0 = result.E_nodes,)) # this is agnostic to the solver
  end
  return settings, result
end

# Global solver.
function solve_full_model_global(settings; impose_E_monotonicity_constraints = true)
  ranges = map(i->(settings.transition_lb[i], settings.transition_ub[i]), 1:length(settings.transition_x0))
  ssr(residuals) = sum(residuals .* residuals)
  result = bboptimize(x -> ssr(residuals_given_E_nodes_interior(impose_E_monotonicity_constraints ? sort(x) : x, settings));
                      SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = settings.transition_iterations)
  return (solution = solve_model_from_E_nodes(best_candidate(result), settings; detailed_solution = true), E_nodes = best_candidate(result))
end

function solve_full_model(settings; impose_E_monotonicity_constraints = true)
   # constraint for increasing E nodes
  function constraints_increasing_E!(h, x, jacobian_t)
    M = length(x)
    # A is a matrix whose ith row has 1 in ith col and -1 in (i+1)th col
    # so ith element of A*x imposes x[i] <= x[i+1]
    # note that this imposes x[end] <= 0 as well as the last row is [0; 0; ...; 0; 1].
    A = LinearAlgebra.Tridiagonal(zeros(M-1),ones(M),-ones(M-1))
    if length(jacobian_t) > 0 # transpose of the Jacobian matrix
        jacobian_t[:] = A'
    end
    h[:] = A*x
  end
   result = solve_system(x -> residuals_given_E_nodes_interior(x, settings), settings.transition_x0;
                    lb = nothing, ub = fill(0.0, length(settings.transition_x0)),
                    constraints_fg! = impose_E_monotonicity_constraints ? constraints_increasing_E! : nothing,
                    algorithm = impose_E_monotonicity_constraints ? :LD_SLSQP : :LD_LBFGS,
                    iterations = settings.transition_iterations)
  return (solution = solve_model_from_E_nodes(result, settings; detailed_solution = true),
          E_nodes = result)
end

function solve_full_model_dfols(settings; user_params = nothing)
  if (settings.transition_iterations == 0)
    result = (x = settings.transition_x0,) # so that result.x is always defined
  else
    result = DFOLS.solve(x -> residuals_given_E_nodes_interior(x, settings), settings.transition_x0, maxfun = settings.transition_iterations, user_params = user_params)
  end
  return (solution = solve_model_from_E_nodes(result.x, settings; detailed_solution = true), E_nodes = result.x, solobj = result)
end

#=
  Auxiliary functions. Organized from most fundamental to least fundamental.
=#

# Returns model solution (i.e., solve_dynamics output) given a set of E nodes. Used by all solvers.
function solve_model_from_E_nodes(E_nodes_interior, settings; detailed_solution = false, interp = LinearInterpolation)
  @unpack T, params_T, stationary_sol_T, Ω_0, E_node_count, entry_residuals_nodes_count, iterations = settings
  δ = params_T.δ
  Ω_T = stationary_sol_T.Ω
  # fix the point at T to be zero and sort candidate if needed
  E_nodes = [-1.0; E_nodes_interior; 0.0] # fix the end points to remove indeterminacy problems
  ts = range(0.0, stop=T, length=length(E_nodes))
  E_hat_interpolation = interp(ts, E_nodes) # might worth trying cubic spline
  E_hat(t) = E_hat_interpolation(t)
  E_hat_integral = quadgk(E_hat, 0, T)[1]
  # Formulate and solve ODEProblem
  M = log(Ω_T/Ω_0) / E_hat_integral # when Ω_T = Ω_0, then M = 0 so that E(t) is constant with δ as expected
  Ω_derivative(Ω,p,t) = M*E_hat(t)*Ω
  Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15) # if this fails, error will be thrown
  Ω(t) = t <= T ? Ω_solution(t) : Ω_solution(T)
  E(t) = M*E_hat(t) + δ
  # solve the dynamics and get the resulting entry_residual vector; if solution is not valid, return Inf
  return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; detailed_solution = detailed_solution)
end

# Call the above and then get the residual. Entry point for DFOLS solver.
function residuals_given_E_nodes_interior(E_nodes_interior, settings)
  solved = try solve_model_from_E_nodes(E_nodes_interior, settings).results catch; return fill(10e20, settings.entry_residuals_nodes_count) end
  entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
  entry_residuals_nodes = range(0, stop = solved.t[end], length = settings.entry_residuals_nodes_count + 2)
  return entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
end