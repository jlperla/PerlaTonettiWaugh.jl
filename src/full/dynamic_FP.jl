# Called at end of each time-step in the DAE solver (mainly to calculate L_tilde for derivatives)
function store_static_callback_FP(u, t, integrator, p)
    @unpack results, parameters, settings, L_1, L_2, ω, Ω, Ξ₁ = p
    P = length(u) - 3
    v_1_t = u[1]
    g_t = u[P+1]
    z_hat_t = u[P+2]
    E_t = u[P+3]
    Ω_t = Ω(t)
    S_t = S(g_t, parameters)
    L_tilde_t = L_tilde(g_t, z_hat_t, Ω_t, E_t, S_t, parameters)
    entry_residual_t = entry_residual(v_1_t, Ξ₁, parameters)
    push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_1 = v_1_t, L_tilde = L_tilde_t, entry_residual = entry_residual_t))
end

function f_FP!(residual,du,u,p,t)
    @unpack results, parameters, settings, L_1, L_2, ω, Ω, Ξ₁ = p
    @unpack z_ex, T = settings
    z = z_ex[2:end-1]
    @unpack ρ, δ, σ, μ, υ, κ, d, ζ = parameters
    residual .= 0
    P = length(residual) - 3 # note u[1:P]  = v(t) in the solution iterators
    v_1 = u[1]
    g = u[P+1]
    z_hat = u[P+2]
    E = u[P+3]

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S, L_tilde, z_bar, π_min, π = static_equilibrium(Ξ₁, v_1, g, z_hat, E, Ω(t), z, parameters)
    L_tilde_log_derivative = 0.0 # Default to Float literal.
    if (t < T)
        t_forward = results[:t][end]
        L_tilde_forward = results[:L_tilde][end]
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

function DAE_outofdomain(u, p, t)
    v_1 = u[1]
    g = u[end - 2]
    z_hat = u[end - 1]
    E = u[end]
    return (v_1 < 0 || g < 0 || z_hat < 1 || E < 0) ? true : false
end

# Calculate the transition dynamics given a fixed Ω(t) function
function solve_dynamics_FP(parameters, stationary_sol_T, settings, Ω)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z_ex, T, tstops = settings
    z = z_ex[2:end-1]
    P = length(z)
    @assert γ ≈ 1 # These are the only supported parameters for the transition dynamics at this point
    @assert η == 0
    T >= 50 ? nothing : @warn "T < 50 can give erratic behavior."

    # unpack the stationary solution
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    L_tilde_T = stationary_sol_T.L_tilde
    Ω_T = stationary_sol_T.Ω

    # Define the results data frame we'll be using and push the stationary onto it.
    results = DataFrame(t = T, g = g_T, z_hat = z_hat_T, Ω = Ω_T, E = δ, v_1 = v_T[1], L_tilde = L_tilde_T)

    # Define intermediate quantitities.
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)

    # initial conditions and parameters for the solver
    u0 = [v_T; g_T; z_hat_T; δ] # i.e. E(T) = δ
    du0 = zeros(P+3)

    p = (results = results, parameters = parameters, settings = settings,
        L_1 = L_1, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, Ω = Ω)

    dae_prob = DAEProblem(f_FP!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false; false])

    # Solve that DAE passing in the callbacks to store "future" L_tilde
    cb = FunctionCallingCallback((u, t, integrator) -> store_static_callback_FP(u, t, integrator, p), tdir = -1, func_start = false) # Callback object.
    sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops, isoutofdomain = DAE_outofdomain)

    # postprocess and return
    sort!(results) # sort results to be forward in time

    return (results = results, # dataframe of results
            sol = sol, # DAE solver output
            dae_parameters = p) # parameters bundle we passed to the DAE
end

# takes an E_hat guess and spits out a new E_hat guess
function dynamics_fixedpoint(E_nodes_interior, params_T, stationary_sol_0, stationary_sol_T, settings)
    # calculate the Ω(t) from the E_hat. Throws away the E (which the DAE calculates)
    @unpack Ω = Ω_E_funcs(E_nodes_interior, stationary_sol_T.Ω, stationary_sol_0.Ω, params_T, settings)
    @unpack δ = params_T
    solution = solve_dynamics_FP(params_T, stationary_sol_T, settings, Ω)

    # turn the DAE solution E into an E_hat (i.e., scale it to [-1, 0])
    filter!(row -> row.t ∈ settings.tstops, solution.results) # cut down to tstops
    E = solution.results[:E]
    E_hat = (E .- δ)/(δ - E[1])

    return E_hat[2:end-1]
end


function solve_transition_FP(parameters, settings)
    @unpack fp_x0, T, tstops = settings
    @unpack d_0, d_T = parameters
    @assert d_0 !== d_T "Will lead to a divide by 0 error."

    # calculate the steady states at 0 and T
    params_T = merge(parameters, (d = d_T,))
    params_0 = merge(parameters, (d = d_0,))
    stationary_sol_T = stationary_numerical(params_T, settings)
    stationary_sol_0 = stationary_numerical(params_0, settings)

    @assert tstops[1] ≈ 0.0
    @assert tstops[end] ≈ T

    fp_sol = fixedpoint(E_nodes_interior -> dynamics_fixedpoint(E_nodes_interior, params_T, stationary_sol_0, stationary_sol_T, settings), fp_x0[2:end-1], iterations = settings.fp_iterations, m = settings.m)
    converged(fp_sol) || @warn "fixed point didn't converge."
    E_nodes_interior = fp_sol.zero

    # regenerate full solution from the E_nodes_interior 
    @unpack Ω = Ω_E_funcs(E_nodes_interior, stationary_sol_T.Ω, stationary_sol_0.Ω, params_T, settings)
    @unpack δ = params_T
    return (dynamics_output = solve_dynamics_FP(params_T, stationary_sol_T, settings, Ω), nodes = E_nodes_interior)
end
