#=
  Functions to solve the model (i.e., return equilibrium quantities and a set of E nodes to interpolate) on a transition from one steady state to another.

  The first set of functions are caller functions, which users would interact with. Then come the auxiliary functions used by the callers.
=#

# Global solver.
function solve_full_model_global(settings; impose_E_monotonicity_constraints = true, front_nodes_appended = nothing)
  # some exception handling on the inputs
  length(settings.transition_lb) == length(settings.transition_ub) || @warn "length(settings.transition_lb) != length(settings.transition_ub)"
  length(settings.transition_x0) == length(settings.transition_ub) || @warn "transition_x0 and transition_ub sizes differ; setting ub to zeros(size(transition_x0))"; settings = merge(settings, (transition_ub = fill(0.0, length(settings.transition_x0)),))
  length(settings.transition_x0) == length(settings.transition_lb) || @warn "transition_x0 and transition_lb sizes differ; setting lb to -1*ones(size(transition_x0))"; settings = merge(settings, (transition_lb = fill(-1.0, length(settings.transition_x0)),))
  length(settings.transition_x0) == length(settings.weights) || @warn "transition_x0 and weights sizes differ; setting weights to default function"; settings = merge(settings, (weights = (x -> x < 10 ? 10. - x : 1.).(1:1:length(settings.transition_x0)),))
  # setup range for bboptimize
  ranges = map(i->(settings.transition_lb[i], settings.transition_ub[i]), 1:length(settings.transition_ub))
  # run solver and process result
  ssr(residuals) = sum(residuals.^2) + settings.transition_penalty_coefficient * sum(max.(-diff(residuals), 0.0)) # add penalty for non-increasing sol
  result = bboptimize(x -> ssr(weighted_residuals_given_E_nodes_interior(impose_E_monotonicity_constraints ? sort(x) : x, settings; front_nodes_appended = front_nodes_appended)); SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = settings.transition_iterations)
  E_nodes_found = best_candidate(result)
  E_nodes_found = impose_E_monotonicity_constraints ? sort(E_nodes_found) : E_nodes_found
  # return
  return (solution = solve_model_from_E_nodes(E_nodes_found, settings; detailed_solution = true), E_nodes = E_nodes_found)
end

function solve_full_model(settings; impose_E_monotonicity_constraints = true, write_csv = false, csvpath = nothing, run_global = true, front_nodes_appended = nothing)
  # preprocessing before NLopt/local solver
  if (run_global && settings.transition_iterations > 0) # run global if required
    result = solve_full_model_global(settings; impose_E_monotonicity_constraints = impose_E_monotonicity_constraints, front_nodes_appended = front_nodes_appended)
    E_nodes = result.E_nodes;
    settings = merge(settings, (transition_x0 = E_nodes, ));
  elseif transition_iterations < 1 # return immediately with the initial condition if that's required
    E_nodes_interior = settings.transition_x0
    return (solution = solve_model_from_E_nodes(E_nodes_interior, settings; detailed_solution = true), E_nodes = E_nodes_interior)
  end
   # linear constraint for increasing E nodes
  function constraints_increasing_E!(h, x, jacobian_t)
    L = length(x)
    # A is a matrix whose ith row has 1 in ith col and -1 in (i+1)th col
    # so ith element of A*x imposes x[i] <= x[i+1]
    # note that this imposes x[end] <= 0 as well as the last row is [0; 0; ...; 0; 1].
    A = LinearAlgebra.Tridiagonal(zeros(L-1),ones(L),-ones(L-1))
    if length(jacobian_t) > 0
      jacobian_t[:] = A' # transpose of the Jacobian matrix
    end
    h[:] = A*x
  end
  # run the local (NLopt) solver
  E_nodes = solve_system(x -> weighted_residuals_given_E_nodes_interior(x, settings; front_nodes_appended = front_nodes_appended), settings.transition_x0;
                    lb = nothing, ub = fill(0.0, length(settings.transition_x0)),
                    constraints_fg! = impose_E_monotonicity_constraints ? constraints_increasing_E! : nothing,
                    algorithm = impose_E_monotonicity_constraints ? :LD_SLSQP : :LD_LBFGS,
                    iterations = settings.transition_iterations)
  # append front_nodes_appended in front if needed
    E_nodes = front_nodes_appended == nothing ? E_nodes : [front_nodes_appended; E_nodes]
    solution = solve_model_from_E_nodes(E_nodes, settings; detailed_solution = true)
  # output caching
    write_csv && CSV.write(csvpath, solution.results)
  # return
  return (solution = solution, E_nodes = E_nodes)
end

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

function weighted_residuals_given_E_nodes_interior(E_nodes_interior, settings; front_nodes_appended = nothing)
  entry_residuals_nodes_count = length(E_nodes_interior)
  # append front_nodes_appended in front if needed
  E_nodes_interior = front_nodes_appended == nothing ? E_nodes_interior : [front_nodes_appended; E_nodes_interior]
  # solve the model using E_nodes_interior received
  solved = try solve_model_from_E_nodes(E_nodes_interior, settings).results catch; return fill(10e20, entry_residuals_nodes_count) end
  entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
  entry_residuals_nodes = range(0, stop = solved.t[end], length = entry_residuals_nodes_count + 2)
  # return weighted residuals
  return settings.weights .* entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
end
