# Continuation solver. With default parameters we get: [-1.1965566608914664,-0.7506441878495663,-0.6160629809004126,-0.4490604063250066,-0.3457782022324279,-0.26611899012195817,-0.1703597837124977,-0.12538756277279978,-0.10765936098467942,-0.057163395444530966,-0.05297636889833756,-0.032260206198113685,-0.028143386570548785,-0.04434566417370368]
function solve_continuation(d_0, d_T; step = 0.005, params = parameter_defaults(), settings = settings_defaults())
  params_0 = merge(params, (d = d_T,)) # parameters to be used at t = 0
  params_T = merge(params, (d = d_T,)) # parameters to be used at t = T
  z_grid = settings.z
  stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
  stationary_sol_T = stationary_numerical(params_T, z_grid) # solution at t = T
  Ω_0 = stationary_sol_0.Ω;
  settings = merge(settings, (params_T = params_T, stationary_sol_T = stationary_sol_T, Ω_0 = Ω_0));
  tempd_0 = params_0.d
  result = 0 # to be overwritten later
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

function solve_full_model_global(settings)
  settings = merge(settings, (iterations = settings.global_transition_iterations,
                              weights = settings.global_transition_weights,
                              sort_candidate = true))
  ranges = map(i->(settings.global_transition_lb[i], settings.global_transition_ub[i]),
              1:length(settings.global_transition_x0))
  result = bboptimize(x -> ssr_given_E_nodes(x, settings); SearchRange = ranges,
                      NumDimensions = length(ranges), MaxSteps = settings.iterations)
  return (solution = solve_with_E_nodes(best_candidate(result), settings; detailed_solution = true),
          E_nodes = best_candidate(result))
end

function solve_full_model_nlopt(settings)
  settings = merge(settings, (sort_candidate = false,))

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

  result = find_zero(x -> residuals_given_E_nodes(x, settings), settings.global_transition_x0; 
                    lb = nothing, ub = fill(0.0, length(settings.global_transition_x0)),
                    constraints_fg! = constraints_increasing_E!, autodiff = :finite)
  return (solution = solve_with_E_nodes(result, settings; detailed_solution = true), 
          E_nodes = result)
end

function solve_full_model_python(settings; user_params = nothing)
  settings = merge(settings, (sort_candidate = false,))
  result = DFOLS.solve(x -> residuals_given_E_nodes(x, settings), settings.global_transition_x0, user_params = user_params)
  return (solution = solve_with_E_nodes(result.x, settings; detailed_solution = true), E_nodes_and_T = result.x, solobj = result)
end

function solve_with_E_nodes(E_nodes, settings; detailed_solution = false, interp = CubicSplineInterpolation)
  @unpack T, params_T, stationary_sol_T, Ω_0, E_node_count, entry_residuals_nodes_count, iterations, sort_candidate = settings

  δ = params_T.δ
  Ω_T = stationary_sol_T.Ω

  # fix the point at T to be zero and sort candidate if needed
  E_nodes = sort_candidate ? [sort(E_nodes); 0.0] : [E_nodes; 0.0]

  # construct Ω and E_nodes
  E_hat_vec_range = E_nodes[end] - E_nodes[1]
  E_hat_vec_scaled = (E_nodes .- E_nodes[1]) ./ E_hat_vec_range .- 1.0
  ts = range(0.0, stop=T, length=length(E_nodes))
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

function residuals_given_E_nodes(E_nodes, settings)
  # solve the dynamics; if solution is not valid, return Inf
  solved = try solve_with_E_nodes(E_nodes, settings).results catch; return fill(10e20, settings.entry_residuals_nodes_count) end
  # get the resulting entry_residual vector
  return residuals_given_solution(solved, settings.entry_residuals_nodes_count)
end

function ssr_given_E_nodes(E_nodes, settings)
  residuals = residuals_given_E_nodes(E_nodes, settings)
  ssr_rooted = sqrt(sum(residuals .* settings.weights .* residuals))
  return ssr_rooted +
          ((settings.global_transition_penalty_coefficient > 0.) ? # add a penalty function for constraints on increasing E
          (settings.global_transition_penalty_coefficient * sum((max.(0.0, diff(E_nodes))).^2)) : 
          0.) 
end
