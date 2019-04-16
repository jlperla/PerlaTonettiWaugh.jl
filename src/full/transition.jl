# Functions to solve the model (i.e., return equilibrium quantities and a set of E nodes to interpolate) on a transition from one steady state to another.
# Global solver.
function solve_transition_global(parameters, settings)
    @assert size(settings.transition_lb) == size(settings.transition_ub) "Lower and upper bound vector lengths disagree."
    @assert size(settings.transition_x0) == size(settings.weights) "x0 and weight vector lengths disagree."
    @assert size(settings.transition_x0) == size(settings.transition_lb) "x0 and bound vector lengths disagree."
    ranges = map(i->(settings.transition_lb[i], settings.transition_ub[i]), 1:length(settings.transition_ub))
    
    # define objective function 
    ssr(residuals) = sum(residuals.^2) + settings.transition_penalty_coefficient * sum(max.(-diff(residuals), 0.0)) # add penalty for non-increasing sol
    bb_objective(x) = (settings.impose_E_monotonicity_constraints) ? ssr(weighted_residuals(sort(x), parameters, settings)) : ssr(weighted_residuals(x, parameters, settings))
    
    # run solver, postprocess, and return 
    result = bboptimize(bb_objective, SearchRange = ranges, NumDimensions = length(ranges), MaxSteps = settings.transition_iterations)
    E_nodes_found = settings.impose_E_monotonicity_constraints ? sort(best_candidate(result)) : best_candidate(result)
    return E_nodes_found
end

function solve_transition(parameters, settings)
    # run global if requested 
    if settings.run_global 
        x0 = solve_transition_global(parameters, settings)
    else 
        x0 = settings.transition_x0 
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
        h .= A*x
    end
    
    # run the local (NLopt) solver
    nlopt_objective(x) =  weighted_residuals(x, parameters, settings)
    E_nodes = solve_system(nlopt_objective, x0;
                    lb = nothing, 
                    ub = fill(0.0, length(x0)),
                    constraints_fg! = settings.impose_E_monotonicity_constraints ? constraints_increasing_E! : nothing,
                    algorithm = settings.impose_E_monotonicity_constraints ? :LD_SLSQP : :LD_LBFGS,
                    iterations = settings.transition_iterations)
    
    # solve, cache, and return 
    solution = solve_model_from_E_nodes(E_nodes, parameters, settings; detailed_solution = true)
    return (solution = solution, results = solution.results, E_nodes = E_nodes)
end

function weighted_residuals(E_nodes_interior, parameters, settings)
    entry_residuals_nodes_count = length(E_nodes_interior)
    
    # solve the model using E_nodes
    solved = try solve_model_from_E_nodes(E_nodes_interior, parameters, settings).results catch; return fill(10e20, entry_residuals_nodes_count) end
    entry_residual_interpolated = LinearInterpolation(solved.t, solved.entry_residual)
    entry_residuals_nodes = range(0, stop = solved.t[end], length = entry_residuals_nodes_count + 2)
    
    return settings.weights .* entry_residual_interpolated.(entry_residuals_nodes[2:(end-1)])
  end

# Returns model solution (i.e., solve_dynamics output) given a set of E nodes. Used by all solvers.
function solve_model_from_E_nodes(E_nodes_interior, parameters, settings; detailed_solution = false)
    @unpack T, interp, z_ex = settings 
    @unpack d_0, d_T = parameters 
    
    # calculate stationary objects
    params_T = merge(parameters, (d = d_T,)) 
    params_0 = merge(parameters, (d = d_0,))
    stationary_sol_T = stationary_numerical(params_T, settings)
    Ω_0 = stationary_numerical(params_0, settings).Ω
    Ω_T = stationary_sol_T.Ω
    
    # formulate and solve the DAEproblem
    @unpack Ω, E = Ω_E_funcs(E_nodes_interior, Ω_T, Ω_0, parameters, settings)
    return solve_dynamics(params_T, stationary_sol_T, settings, Ω, E; detailed_solution = detailed_solution)
end

function Ω_E_funcs(E_nodes_interior, Ω_T, Ω_0, parameters, settings)
    @unpack δ = parameters
    @unpack T, interp = settings
    E_nodes = [-1.0; E_nodes_interior; 0.0] # See footnote 19.  Without pinning a corner there is an extra degree of freedom in the scale
    ts = range(0.0, T, length=length(E_nodes))

    # calculate E 
    E_hat_interpolation = interp(ts, E_nodes) # might worth trying cubic spline
    E_hat(t) = E_hat_interpolation(t)
    E_hat_integral = quadgk(E_hat, 0, T)[1] # (42)
    Q = log(Ω_T/Ω_0) / E_hat_integral # (42) when Ω_T = Ω_0, then Q = 0 so that E(t) is constant with δ as expected
    E(t) = Q*E_hat(t) + δ # (43)

    # calculate Ω
    Ω_derivative(Ω,p,t) = Q*E_hat(t)*Ω # (44)
    Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative, Ω_0, (0.0, T)), reltol = 1e-15) # if this fails, error will be thrown
    Ω(t) = t <= T ? Ω_solution(t) : Ω_solution(T)

    return (Ω = Ω, E = E)
end 
