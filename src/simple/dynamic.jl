# DAE kernel, defined for each time ∈ t.
function f!_simple(resid,du,u,p,t)
    @unpack Ξ₁, L_1, L_2, z, r, μ, υ, π, T, P, x, ω = p

    v_t = u[1:P]
    g_t = u[P+1]
    A = (r(t) - μ - υ^2/2)*I - (μ + υ^2 - g_t)*L_1 - υ^2/2 * L_2 # (19)

    resid[1:P] .= A * v_t - π.(t, z) # (28)
    resid[1:P] .-= du[1:P] # discretized system of ODE for v (28)
    resid[P+1] = Ξ₁*v_t[1] - dot(ω, v_t) + x(t) # value matching condition (29)
end

# DAE constructor
function simpleDAE(parameters, settings)
    @unpack μ, υ, θ, r, x, π = parameters
    @unpack z_ex, T = settings

    @assert z_ex[1] == 0.0 && issorted(z_ex) # validate grid
    z = z_ex[2:end-1]
    P = length(z)

    ω = ω_weights(z_ex, θ, 1)

    bc = (Mixed(1), Mixed(1)) # boundary conditions for differential operators
    L_1_minus = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)

    stationary = stationary_numerical_simple(parameters, settings)
    g_T = stationary.g
    v_T = stationary.v

    Ξ₁ = 1/(1 - 1*(z[1] - 0.0)) # (24)
    p = (Ξ₁ = Ξ₁, L_1 = L_1_minus, L_2 = L_2, z = z, r = r, υ = υ, π = π, T = T, μ = μ, g_T = g_T, P = P, x = x, ω = ω)

    u_T = [v_T; g_T]
    resid_M1 = zeros(P+1)
    return DAEProblem(f!_simple, resid_M1, u_T, (T, 0.0), differential_vars = [fill(true, P); false], p)
end

function solve_simple_transition(parameters, settings)
    ode_prob = simpleDAE(parameters, settings)
    @unpack ω, x, Ξ₁, P = ode_prob.p
    @unpack ts, ode_solve_algorithm = settings
    sol = Sundials.solve(ode_prob, ode_solve_algorithm, tstops = ts)

    # Calculate the residual at each time point
    residuals = zeros(length(ts))
    v_ts = zeros(P, length(ts))
    g_ts = zeros(length(ts))
    for (i, t) in enumerate(ts)
        v_t = sol(t)[1:P] # i.e., the value function at the point-in-time.
        residuals[i] = Ξ₁*v_t[1] + x(t) - dot(ω, v_t) # value matching condition (29)
        v_ts[:,i] = v_t
        g_ts[i] = sol(t)[end]
    end
    return (residuals = residuals, v_ts = v_ts, g_ts = g_ts)
end
