# Called at end of each time-step in the DAE solver (mainly to calculate L_tilde for derivatives)
function store_static_callback(u, t, integrator, p)
    @unpack results, parameters, settings, L_1, L_2, ω, E, Ω, Ξ₁ = p
    P = length(u) - 2
    g_t = u[P+1]
    z_hat_t = u[P+2]
    Ω_t = Ω(t)
    E_t = E(t)
    v_1_t = u[1]
    S_t = S(g_t, parameters)
    L_tilde_t = L_tilde(g_t, z_hat_t, Ω_t, E_t, S_t, parameters)
    push!(results, (t = t, g = g_t, z_hat = z_hat_t, Ω = Ω_t, E = E_t, v_1 = v_1_t, L_tilde = L_tilde_t))
end

function f!(residual,du,u,p,t)
    @unpack results, parameters, settings, L_1, L_2, ω, E, Ω, Ξ₁ = p
    @unpack z, T = settings
    @unpack ρ, δ, σ, μ, υ, κ, d, ζ = parameters
    residual .= 0
    P = length(residual) - 2 # note u[1:P]  = v(t) in the solution iterators
    g = u[P+1]
    z_hat = u[P+2]

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S, L_tilde, z_bar, π_min, π = static_equilibrium(Ξ₁, u[1], g, z_hat, E(t), Ω(t), z, parameters)
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
    residual[P+1] = Ξ₁*u[1] + ζ - dot(ω, u[1:P])  # (54)
    residual[P+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min  # (55)
end

# Calculate the transition dynamics given a fixed Ω(t), E(t) function
# To calculate the full equilibrium, you can solve until Ω(t), E(t) fulfill free-entry conditions
function solve_dynamics(parameters, stationary_sol_T, settings, Ω, E)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z, z_ex, T, tstops = settings
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
    P = length(z)
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)

    # initial conditions and parameters for the solver
    u0 = [v_T; g_T; z_hat_T]
    du0 = zeros(P+2)

    p = (results = results, parameters = parameters, settings = settings,
        L_1 = L_1, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, E = E, Ω = Ω)

    dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false])

    # Solve that DAE passing in the callbacks to store "future" L_tilde
    cb = FunctionCallingCallback((u, t, integrator) -> store_static_callback(u, t, integrator, p), tdir = -1, func_start = false) # Callback object.
    if (tstops == nothing)
        sol = DifferentialEquations.solve(dae_prob, callback = cb)
    else
        sol = DifferentialEquations.solve(dae_prob, callback = cb, tstops = tstops)
    end

    # postprocess and return
    sort!(results) # sort results to be forward in time
    results = @transform(results, entry_residual = entry_residual.(:v_1, Ref(Ξ₁), Ref(parameters)))

    return (results = results, # dataframe of results
            sol = sol, # DAE solver output
            dae_parameters = p) # parameters bundle we passed to the DAE
end
