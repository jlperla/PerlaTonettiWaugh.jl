# Called at end of each time-step in the DAE solver (mainly to calculate L_tilde for derivatives)
function store_static_callback(u, t, integrator, p)
    @unpack results, parameters, settings, L_1, L_2, ω, Ω, Ξ₁ = p
    P = length(u) - 3
    v_1_t = u[1]
    g_t = u[P+1]
    z_hat_t = max(u[P+2], 1.)
    E_t = u[P+3]
    Ω_t = Ω(t)
    S_t = S(g_t, parameters)
    L_tilde_t = L_tilde(g_t, z_hat_t, Ω_t, E_t, S_t, parameters)
    entry_residual_t = entry_residual(v_1_t, Ξ₁, parameters)
    push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_1 = v_1_t, L_tilde = L_tilde_t, entry_residual = entry_residual_t))
end

# Dynamic kernel
function f!(residual,du,u,p,t)
    @unpack results, parameters, settings, L_1, L_2, ω, Ω, Ξ₁, z = p
    @unpack z_ex, T = settings
    @unpack ρ, δ, σ, μ, υ, κ, d, ζ = parameters
    residual .= 0
    P = length(residual) - 3 # note u[1:P]  = v(t) in the solution iterators
    v_1 = u[1]
    g = u[P+1]
    z_hat = max(u[P+2], 1.)
    E = u[P+3]

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S, L_tilde, z_bar, π_min, π = static_equilibrium(Ξ₁, v_1, g, z_hat, E, Ω(t), z, parameters)
    L_tilde_log_derivative = 0.0 # Default to Float literal.
    if (t < T)
        t_forward = results.t[end]
        L_tilde_forward = results.L_tilde[end]
        L_tilde_log_derivative = (log(1 - L_tilde_forward) - log(1 - L_tilde))/(t_forward - t) # See note under (40)
    end

    # Combine (40) with (52) to get A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2
    # Then building the A(t) * v(t) - v'(t) residual directly for the non-algebraic equations
    residual[1:P] = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:P] # (40) into (52) then (53)
    residual[1:P] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:P] # (52) 2nd term in (53)
    residual[1:P] .-= (υ^2/2)*L_2*u[1:P] # (52) third term in (53)
    residual[1:P] .-= π # (53) final term
    residual[1:P] .-= du[1:P] # (53) subtracting the v'(t) to form the residual
    residual[P+1] = Ξ₁*v_1 + ζ - dot(ω, u[1:P])  # (54)
    residual[P+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min  # (55)
    residual[P+3] = entry_residual(v_1, Ξ₁, parameters)  # (56)
end

# Calculate the transition dynamics given a fixed Ω(t) function
function solve_dynamics(parameters, stationary_sol_T, settings, Ω)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z_ex, T, tstops = settings
    z = z_ex[2:end-1]
    P = length(z)
    @assert γ ≈ 1 # These are the only supported parameters for the transition dynamics at this point
    @assert η == 0

    # unpack the stationary solution
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    L_tilde_T = stationary_sol_T.L_tilde
    Ω_T = stationary_sol_T.Ω
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)
    entry_residual_T = entry_residual(v_T[1], Ξ₁, parameters)

    # Define the results data frame we'll be using and push the stationary onto it.
    results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_1 = v_T[1], L_tilde = L_tilde_T, entry_residual = entry_residual_T)

    # Define intermediate quantitities.
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(ξ = σ-1), Mixed(ξ = σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)

    # initial conditions and parameters for the solver (under backwards time)
    u0 = [v_T; g_T; z_hat_T; δ] # i.e. E(T) = δ
    du0 = zeros(P+3)

    p = (results = results, parameters = parameters, settings = settings,
        L_1 = L_1, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, Ω = Ω, z = z)

    dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false; false])

    # Solve that DAE passing in the callbacks to store "future" L_tilde
    cb = FunctionCallingCallback((u, t, integrator) -> store_static_callback(u, t, integrator, p), tdir = -1, func_start = false) # Callback object.
    sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops)

    # postprocess and return
    sort!(results) # sort results to be forward in time

    return (results = results, # dataframe of results
            sol = sol, # DAE solver output
            dae_parameters = p) # parameters bundle we passed to the DAE
end

# takes an E_hat guess and spits out a new E_hat guess
function dynamics_fixedpoint(E_nodes_interior, params_T, stationary_sol_0, stationary_sol_T, settings)
    # calculate the Ω(t) from the E_hat. Throws away the E (which the DAE calculates)
    Ω = Ω_from_E_hat(E_nodes_interior, stationary_sol_T.Ω, stationary_sol_0.Ω, params_T, settings)
    @unpack δ = params_T
    solution = solve_dynamics(params_T, stationary_sol_T, settings, Ω)

    # turn the DAE solution E into an E_hat (i.e., scale it to [-1, 0])
    filter!(row -> row.t ∈ settings.tstops, solution.results) # cut down to tstops
    E = solution.results[:E]
    E_hat = (E .- δ)/(δ - E[1]) # (43)

    return E_hat[2:end-1]
end

function solve_transition(parameters, settings)
    @unpack T, tstops, fixedpoint_beta, fixedpoint_ftol, fixedpoint_iterations, fixedpoint_show_trace, fixedpoint_m, fixedpoint_x0  = settings
    x0 = fixedpoint_x0(parameters, settings)
    @unpack d_0, d_T = parameters
    @assert d_0 !== d_T "Will lead to a divide by 0 error."

    # calculate the steady states at 0 and T
    params_T = merge(parameters, (d = d_T,))
    params_0 = merge(parameters, (d = d_0,))
    stationary_sol_T = stationary_numerical(params_T, settings)
    stationary_sol_0 = stationary_numerical(params_0, settings)

    @assert tstops[1] ≈ 0.0
    @assert tstops[end] ≈ T

    fp_sol = fixedpoint(E_nodes_interior -> dynamics_fixedpoint(E_nodes_interior, params_T, stationary_sol_0, stationary_sol_T, settings),
                        x0,
                        iterations = fixedpoint_iterations,
                        ftol = fixedpoint_ftol,
                        m = fixedpoint_m,
                        beta = fixedpoint_beta,
                        show_trace = fixedpoint_show_trace)
    converged(fp_sol) || @warn "fixed point didn't converge."
    E_nodes_interior = fp_sol.zero

    # regenerate full solution from the E_nodes_interior
    Ω = Ω_from_E_hat(E_nodes_interior, stationary_sol_T.Ω, stationary_sol_0.Ω, params_T, settings)
    @unpack δ = params_T
    solve_dynamics_output = solve_dynamics(params_T, stationary_sol_T, settings, Ω)
    return (solve_dynamics_output = solve_dynamics_output,
            nodes = E_nodes_interior,
            results = prepare_results(solve_dynamics_output, stationary_sol_T, stationary_sol_0))
end

function prepare_results(solve_dynamics_output, stationary_T, stationary_0)
    @unpack results, sol, dae_parameters = solve_dynamics_output
    @unpack Ω, parameters, settings, Ξ₁, ω = dae_parameters
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z_ex, T, tstops = settings
    z = z_ex[2:end-1]
    P = length(z)

    # time-T values for welfare integrals
    L_tilde_T = stationary_T.L_tilde
    g_T = stationary_T.g
    z_hat_T = stationary_T.z_hat
    Ω_T = stationary_T.Ω
    log_c_T = log(c(L_tilde_T, Ω_T, z_bar(z_hat_T, Ω_T, parameters)))

    # interpolates for welfare integrals
    interpolated_L_tilde(t) = L_tilde(sol(t)[P+1], sol(t)[P+2], Ω(t), sol(t)[P+3], S(sol(t)[P+1], parameters), parameters) # recall that the P+1-th element of sol is g, P+2-nd is z_hat, P+3-rd is E
    interpolated_z_bar(t) = z_bar(sol(t)[P+2], Ω(t), parameters)
    interpolated_c(t) = c(interpolated_L_tilde(t), Ω(t), interpolated_z_bar(t))
    log_c(t) = log(interpolated_c(t))

    # welfare rows for dataframe
    log_M(t) = quadgk(t -> sol(t)[P+1], 0, t)[1]
    U(t) = quadgk(τ -> exp(-ρ*τ)*(log_M(t+τ) + log_c(t+τ)), 0, (T-t))[1] + exp(-ρ*(T-t))*(g_T + ρ*(log_c_T + g_T * T))/(ρ^2) # (C.83)

    # build out dataframe
    results = @transform(results, Ω = Ω.(:t))
    results = @transform(results, λ_ii = λ_ii.(:z_hat, Ref(parameters)))
    results = @transform(results, S = S.(:g, Ref(parameters)))
    results = @transform(results, z_bar = z_bar.(:z_hat, :Ω, Ref(parameters)))
    results = @transform(results, c = c.(:L_tilde, :Ω, :z_bar))
    results = @transform(results, π_min = π_min.(:L_tilde, :z_bar, Ref(parameters)))
    results = @transform(results, log_M = log_M.(:t))
    results = @transform(results, U = U.(:t))
    results = @transform(results, π_rat = π_rat.(:z_hat, Ref(parameters)))
    results = @transform(results, L_tilde_a = L_tilde_a.(:Ω, :S, Ref(parameters)))
    results = @transform(results, L_tilde_x = L_tilde_x.(:z_hat, :Ω, Ref(parameters)))
    results = @transform(results, L_tilde_E = L_tilde_E.(:Ω, :E, Ref(parameters)))
    results = @transform(results, w = w.(:z_bar, Ref(parameters)))

    results.r = ones(Float64, nrow(results)) # filler, to be overwritten
    for i in 1:nrow(results)
        t = results.t[i]
        c = results.c[i]
        g = results.g[i]
        log_c_forward = (i < nrow(results)) ? (log(results.c[i+1]) - log(c))/(results.t[i+1] - t) : 0.0
        if (i != nrow(results))
            @assert results.t[i+1] > t
        end
        results.r[i] = ρ +  δ + γ*(g + log_c_forward) # (C.56)
    end

    # prepend pre-shock steadystate to dataframe
    reduced_stationary = delete(stationary_0, :F, :a, :b, :ν, :x, :y, :U_bar, :v_tilde)
    stationary_data = merge(reduced_stationary, (E = δ, U = stationary_0.U_bar, v_1 = stationary_0.v_tilde[1], entry_residual = Ξ₁*stationary_0.v_tilde[1] - dot(stationary_0.v_tilde, ω) + ζ))
    ts_0 = settings.pre_shock_times
    for t in ts_0
        push!(results, merge(stationary_data, (t = t, log_M = t * stationary_data.g)))
    end

    return sort!(results) # sort results again to keep time in the right order
end

function Ω_from_E_hat(E_nodes_interior, Ω_T, Ω_0, parameters, settings)
    @unpack δ = parameters
    @unpack T, tstops, interp = settings
    E_nodes = [-1.0; E_nodes_interior; 0.0] # See footnote 19.  Without pinning a corner there is an extra degree of freedom in the scale
    ts = range(0.0, T, length=length(E_nodes))

    E_hat_interpolation = interp(ts, E_nodes) # might worth trying cubic spline
    E_hat(t) = E_hat_interpolation(t)
    E_hat_integral = quadgk(E_hat, 0, T)[1] # (42)
    Q = log(Ω_T/Ω_0) / E_hat_integral # (42) when Ω_T = Ω_0, then Q = 0 so that E(t) is constant with δ as expected
    E(t) = Q*E_hat(t) + δ # (43)

    Ω_derivative(Ω,p,t) = Q*E_hat(t)*Ω # (44)
    Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative, Ω_0, (0.0, T)), reltol = 1e-15, tstops = tstops) # if this fails, error will be thrown

    function Ω(t)
        if t > T
            return Ω_T
        elseif t > 0 && t <= T
            return Ω_solution(t)
        else # useful for fixed point approach
            return Ω_0
        end
    end

    return Ω
end
