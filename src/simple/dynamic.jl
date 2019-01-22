# DAE kernel, defined for each time ∈ t.
function f!_simple(resid,du,u,p,t)
    @unpack L_1, L_2, z, r, μ, g, υ, π, T, P, x, ω = p
    # Carry out calculations.
    v_t = u[1:P]
    g_t = u[P+1]
    A = (r(t) - μ - υ^2/2)*I - (μ + υ^2 - g_t)*L_1 - υ^2/2 * L_2 # (17)
    # Calculate residuals.
    resid[1:P] .= A * v_t - π.(t, z) # (22)
    resid[1:P] .-= du[1:P] # discretized system of ODE for v (22)
    resid[P+1] = v_t[1] + x(t) - dot(ω, v_t) # value matching condition (23)
end

# DAE constructor
function simpleDAE(params, settings)
    # Unpack necessary objects.
    @unpack μ, υ, θ, r, x, π = params
    @unpack z, T, g = settings
    P = length(z)
    # Quadrature weighting
    ω = ω_weights(z, θ, 1)
    # Discretize the operator.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, 1) # L_1_minus ≡ L_1 is the only one we use.
    # Calculate the stationary solution.
    r_T = r(T)
    g_T = g(T)
    A_T = (r_T - μ - υ^2/2)*I - (μ + υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (17)
    v_T = A_T \ π.(T, z) # (24)
    # Bundle as before.
    p = (L_1 = L_1_minus, L_2 = L_2, z = z, g = g, r = r, υ = υ, π = π, T = T, μ = μ, g_T = g_T, P = P, x = x, ω = ω)
    # Other objects
    u_T = [v_T; g_T]
    resid_M1 = zeros(P+1)
    return DAEProblem(f!_simple, resid_M1, u_T, (T, 0.0), differential_vars = [fill(true, P); false], p)
end

# Calculate residuals given diffeq problem and primitives
function calculate_residuals(ode_prob, x, ω, ode_solve_algorithm, ts) # To keep the params consistent with other tuples.
    # Solve ode
    sol = Sundials.solve(ode_prob, ode_solve_algorithm, tstops = ts)
    P = length(ω)
    # Calculate the residual at each time point
    residuals = zeros(length(ts))
    v_ts = zeros(P, length(ts))
    g_ts = zeros(length(ts))
    for (i, t) in enumerate(ts)
        v_t = sol(t)[1:P] # i.e., the value function at the point.
        residuals[i] = v_t[1] + x(t) - dot(ω, v_t) # value matching condition (23)
        v_ts[:,i] = v_t
        g_ts[i] = sol(t)[end]
    end
    return (residuals = residuals, v_ts = v_ts, g_ts = g_ts)
end
