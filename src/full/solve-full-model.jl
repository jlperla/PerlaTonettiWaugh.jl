function solve_full_model_global(solution0, settings)
  result = bboptimize(x -> ssr_given_candidate(x, settings); SearchRange = settings.ranges, 
                      NumDimensions = length(settings.ranges), MaxSteps = settings.iterations)
  solution = best_candidate(result)

  return solve_with_candidate(solution, settings; detailed_solution = true)
end

function solve_with_candidate(candidate, settings; detailed_solution = false)
  @unpack params_T, stationary_sol_T, Ω_0, E_node_count, entry_residuals_nodes_count, weights, ranges, iterations, sort_candidate = settings
  δ = params_T.δ    
  Ω_T = stationary_sol_T.Ω

  T = candidate[end]

  # fix the point at T to be zero and sort candidate if needed
  candidate = sort_candidate ? [sort(candidate[1:(end-1)]); 0.0] : [candidate[1:(end-1)]; 0.0] 
  
  # construct Ω and E
  E_hat_vec_range = candidate[end] - candidate[1]
  E_hat_vec_scaled = (candidate .- candidate[1]) ./ E_hat_vec_range .- 1.0 
  ts = range(0.0, stop=T, length=length(candidate))
  E_hat_interpolation = CubicSplineInterpolation(ts, E_hat_vec_scaled) # might worth trying cubic spline
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
  solved = try solve_with_candidate(candidate, settings).results catch; return fill(Inf, settings.entry_residuals_nodes_count) end
  # get the resulting entry_residual vector
  return residuals_given_solution(solved, settings.entry_residuals_nodes_count)
end

function ssr_given_candidate(candidate, settings) 
  residuals = residuals_given_candidate(candidate, settings)
  # solve the dynamics; if solution is not valid, return Inf
  return (sqrt(sum(residuals .* settings.weights .* residuals)))
end
