#=
  Functions to solve the model (i.e., return equilibrium quantities and a set of E nodes to interpolate) on a transition from one steady state to another.

  The first set of functions are caller functions, which users would interact with. Then come the auxiliary functions used by the callers.
=#

# Continuation solver. Essentially traces out a smooth path from some d_T to d_0, using the result from one step as the x0 for the next.
# With default parameters we get: [-1.1965566608914664,-0.7506441878495663,-0.6160629809004126,-0.4490604063250066,-0.3457782022324279,-0.26611899012195817,-0.1703597837124977,-0.12538756277279978,-0.10765936098467942,-0.057163395444530966,-0.05297636889833756,-0.032260206198113685,-0.028143386570548785,-0.04434566417370368]
function solve_continuation(d_0, d_T; step = 0.005, params = parameter_defaults(), settings = settings_defaults())
  params_0 = merge(params, (d = d_T,)) # parameters to be used at t = 0
  params_T = merge(params, (d = d_T,)) # parameters to be used at t = T
  z_grid = settings.z
  stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
  stationary_sol_T = stationary_numerical(params_T, z_grid) # solution at t = T
  Ω_0 = stationary_sol_0.Ω;
  settings = merge(settings, (params_T = params_T, stationary_sol_T = stationary_sol_T, Ω_0 = Ω_0));
  tempd_0 = params_0.d
  result = 0 # to be overwritten by loop
  @softscope while tempd_0 <= d_0
    tempd_0 += step
    params_0 = merge(params, (d = tempd_0,))
    Ω_0 = stationary_numerical(params_0, z_grid).Ω
    settings = merge(settings, (Ω_0 = Ω_0,))
    result = solve_full_model_python(settings)
    settings = merge(settings, (global_transition_x0 = result.solobj.x,))
  end
  return settings, result
end

# Global solver.
function solve_full_model_global(settings)
  settings = merge(settings, (iterations = settings.global_transition_iterations, weights = settings.global_transition_weights, sort_candidate = true))
  ranges = map(i->(settings.global_transition_lb[i], settings.global_transition_ub[i]), 1:length(settings.global_transition_x0))
  result = bboptimize(x -> ssr_given_E_nodes(x, settings); SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = settings.iterations)
  return (solution = solve_model_from_E_nodes(best_candidate(result), settings; detailed_solution = true), E_nodes = best_candidate(result))
end

# DFOLS (Derivative-Free Optimization for Least Squares) solver.
function solve_full_model_python(settings; user_params = nothing)
  settings = merge(settings, (sort_candidate = false,))
  result = DFOLS.solve(x -> residuals_given_E_nodes(x, settings), settings.global_transition_x0, user_params = user_params)
  return (solution = solve_model_from_E_nodes(result.x, settings; detailed_solution = true), E_nodes_and_T = result.x, solobj = result)
end

#=
  Auxiliary functions. Organized from most fundamental to least fundamental.
=#

# Returns model solution (i.e., solve_dynamics output) given a set of E nodes. Used by all solvers.
function solve_model_from_E_nodes(E_nodes, settings; detailed_solution = false, interp = CubicSplineInterpolation)
  @unpack T, params_T, stationary_sol_T, Ω_0, E_node_count, entry_residuals_nodes_count, iterations, sort_candidate = settings
  δ = params_T.δ
  Ω_T = stationary_sol_T.Ω
  # fix the point at T to be zero and sort candidate if needed
  E_nodes = sort_candidate ? [sort(E_nodes); 0.0] : [E_nodes; 0.0]
  # construct Ω and E_nodes
  E_hat_vec_range = E_nodes[end] - E_nodes[1]
  E_hat_vec_scaled = (E_hat_vec_range != 0) ? (E_nodes .- E_nodes[1]) ./ E_hat_vec_range .- 1.0 : zeros(length(E_nodes))
  ts = range(0.0, stop=T, length=length(E_nodes))
  E_hat_interpolation = interp(ts, E_hat_vec_scaled) # might worth trying cubic spline
  E_hat(t) = E_hat_interpolation(t)
  # Formulate and solve ODEProblem
  M = log(Ω_T/Ω_0) / quadgk(E_hat, 0, T)[1]
  Ω_derivative(Ω,p,t) = M*E_hat(t)*Ω
  Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15) # if this fails, error will be thrown
  Ω(t) = t <= T ? Ω_solution(t) : Ω_solution(T)
  E(t) = M*E_hat(t) + δ
  # solve the dynamics and get the resulting entry_residual vector; if solution is not valid, return Inf
  return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; detailed_solution = detailed_solution)
end

# Call the above and then get the residual. Entry point for DFOLS solver.
function residuals_given_E_nodes(E_nodes, settings)
  solved = try solve_model_from_E_nodes(E_nodes, settings).results catch; return fill(10e20, settings.entry_residuals_nodes_count) end
  entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
  entry_residuals_nodes = range(0, stop = solved.t[end], length = settings.entry_residuals_nodes_count + 2)
  return entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
end

# Call the above two functions and get the SSR. Entry point for global solver.
function ssr_given_E_nodes(E_nodes, settings)
  residuals = residuals_given_E_nodes(E_nodes, settings)
  ssr_rooted = sqrt(sum(residuals .* settings.weights .* residuals))
  return ssr_rooted +
          ((settings.global_transition_penalty_coefficient > 0.) ? # add a penalty function for constraints on increasing E
          (settings.global_transition_penalty_coefficient * sum((max.(0.0, diff(E_nodes))).^2)) : 0.)
end
