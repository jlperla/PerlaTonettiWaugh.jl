function solve_full_model_global(settings)
  settings = merge(settings, (iterations = settings.global_transition_iterations, 
                              weights = settings.global_transition_weights,
                              sort_candidate = true))
  ranges = map(i->(settings.global_transition_lb[i], settings.global_transition_ub[i]), 
              1:length(settings.global_transition_x0))
  result = bboptimize(x -> ssr_given_candidate(x, settings); SearchRange = ranges,
                      NumDimensions = length(ranges), MaxSteps = settings.iterations)
  return (solution = solve_with_candidate(best_candidate(result), settings; detailed_solution = true),
          E_nodes_and_T = best_candidate(result))
end

function solve_full_model_python(settings; user_params = nothing)
  settings = merge(settings, (sort_candidate = false,))
  result = DFOLS.solve(x -> residuals_given_candidate(x, settings), settings.global_transition_x0, user_params = user_params)
  return (solution = solve_with_candidate(result.x, settings; detailed_solution = true,
          E_nodes_and_T = result.x))
end

function solve_with_candidate(candidate, settings; detailed_solution = false, interp = CubicSplineInterpolation)
  @unpack T, params_T, stationary_sol_T, Ω_0, E_node_count, entry_residuals_nodes_count, weights, iterations, sort_candidate = settings

  δ = params_T.δ
  Ω_T = stationary_sol_T.Ω

  # fix the point at T to be zero and sort candidate if needed
  candidate = sort_candidate ? [sort(candidate[1:end]); 0.0] : [candidate[1:end]; 0.0]

  # construct Ω and E
  E_hat_vec_range = candidate[end] - candidate[1]
  E_hat_vec_scaled = (candidate .- candidate[1]) ./ E_hat_vec_range .- 1.0
  ts = range(0.0, stop=T, length=length(candidate))
  E_hat_interpolation = interp(ts, E_hat_vec_scaled) # might worth trying cubic spline
  E_hat(t) = E_hat_interpolation(t)

  M = log(Ω_T/Ω_0) / quadgk(E_hat, 0, T)[1]
  Ω_derivative(Ω,p,t) = M*E_hat(t)*Ω
  Ω_solution = try DifferentialEquations.solve(ODEProblem(Ω_derivative,Ω_0,(0.0, T)), reltol = 1e-15) catch; return Inf end
  Ω(t) = t <= T ? Ω_solution(t) : Ω_solution(T)
  E(t) = M*E_hat(t) + δ

  # solve the dynamics and get the resulting entry_residual vector; if solution is not valid, return Inf
  return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; stopwithf! = false, detailed_solution = detailed_solution)
end

function residuals_given_solution(solved, entry_residuals_nodes_count)
  # interpolate on returned entry_residual
  entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)

  # evaluate entry_residual on entry_residual_nodes, return the norm
  entry_residuals_nodes = range(0, stop = solved.t[end], length = entry_residuals_nodes_count + 2)

  # returns the vector of residuals
  return entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
end

function residuals_given_candidate(candidate, settings)
  # solve the dynamics; if solution is not valid, return Inf
  solved = try solve_with_candidate(candidate, settings).results catch; return fill(10e20, settings.entry_residuals_nodes_count) end
  # get the resulting entry_residual vector
  return residuals_given_solution(solved, settings.entry_residuals_nodes_count)
end

function ssr_given_candidate(candidate, settings)
  residuals = residuals_given_candidate(candidate, settings)
  # solve the dynamics; if solution is not valid, return Inf
  return (sqrt(sum(residuals .* settings.weights .* residuals)))
end
