function solve_full_model_global(solution0, params_T, stationary_sol_T, Ω_0, settings; stopwithf! = false, detailed_solution = true)
    @unpack δ = params_T
    @unpack E_node_count, entry_residuals_nodes_count, weights, ranges, iterations = settings
  
    Ω_T = stationary_sol_T.Ω
  
    function solve_with_candidate(candidate)
      candidate = [candidate...] # if candidate is a tuple, convert it to an array
      T = candidate[end]
      
      candidate = [sort(candidate[1:(end-1)]); 0.0] # fix the point at T to be zero
      
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
      return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω, E; stopwithf! = stopwithf!, detailed_solution = detailed_solution)
    end
  
    function residuals_given_solution(solved, entry_residuals_nodes_count)
      # interpolate on returned entry_residual
      entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
    
      # evaluate entry_residual on entry_residual_nodes, return the norm
      entry_residuals_nodes = range(0, stop = solved.t[end], length = entry_residuals_nodes_count + 2)
      
      # returns the vector of residuals
      return entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
    end
    
    function evaluate_candidate(candidate)
      # solve the dynamics; if solution is not valid, return Inf
      solved = try solve_with_candidate(candidate).results catch; return Inf end
      # get the resulting entry_residual vector
      residuals = residuals_given_solution(solved, entry_residuals_nodes_count)
      return (sqrt(sum(residuals .* weights .* residuals)))
    end
  
    result = bboptimize(evaluate_candidate; SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = iterations)
    solution = best_candidate(result)
  
    return solve_with_candidate(solution)
end